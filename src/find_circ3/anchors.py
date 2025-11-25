# src/find_circ3/anchors.py

"""
unmapped2anchors3: Python 3 port of the legacy unmapped2anchors.py

Usage (BAM mode, typical):

    samtools view -hf 4 aln.bam | samtools view -Sb - > unmapped.bam
    find-circ3-anchors unmapped.bam -a 20 -q 5 > anchors.fastq

This module reproduces the original behaviour:

- Takes unmapped reads (or FASTA/FASTQ/reads text) as input.
- Applies minimal anchor-quality filtering.
- Emits A/B anchors in FASTQ, with optional permutation / reverse modes.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable, Iterator, Tuple

import click
import numpy as np
import pysam

# --- Nucleotide complement helpers (ported from legacy) -------------------- #

COMPLEMENT = {
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
    "k": "m",
    "m": "k",
    "r": "y",
    "y": "r",
    "s": "s",
    "w": "w",
    "b": "v",
    "v": "b",
    "h": "d",
    "d": "h",
    "n": "n",
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "K": "M",
    "M": "K",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "B": "V",
    "V": "B",
    "H": "D",
    "D": "H",
    "N": "N",
}


def complement(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in seq)


def rev_comp(seq: str) -> str:
    return complement(seq)[::-1]


# --- FASTA / FASTQ chunkers (ported) --------------------------------------- #


def fasta_chunks(
    lines: Iterable[str], strip: bool = True, fuse: bool = True
) -> Iterator[Tuple[str, str]]:
    """
    Yield (header, sequence) tuples from a FASTA stream.
    """
    chunk = ""
    data = []

    for l in lines:
        if l.startswith("#"):
            continue
        if l.startswith(">"):
            if data and chunk:
                yield chunk, "".join(data)
                if strip:
                    data = []
                else:
                    data.append(l)
            chunk = l[1:].rstrip() if strip else l
        else:
            data.append(l.strip())

    if data and chunk:
        yield chunk, "".join(data)


def qfa_chunks(lines: Iterable[str]):
    """
    Iterate over FASTQ lines, yielding simple objects with
    .name, .seq, .qual (legacy qfa_tuple behaviour).
    """
    from collections import namedtuple

    QFA = namedtuple("qfa_tuple", "name,seq,qual")
    it = iter(lines)

    try:
        while True:
            name = next(it).rstrip()
            if not name.startswith("@"):
                # tolerate blank lines etc.
                continue
            seq = next(it).rstrip()
            plus = next(it).rstrip()
            qual = next(it).rstrip()
            yield QFA(name[1:], seq, qual)  # drop leading "@"
    except StopIteration:
        return


# --- Core anchor emission logic -------------------------------------------- #


def _build_funcs(rev_mode: str):
    """
    Port of legacy funcs dict:

    read_f, A_f, B_f = funcs[rev_mode]
    """
    def passthru(x: str) -> str:
        return x

    def reverse(x: str) -> str:
        return x[::-1]

    funcs = {
        "A": (passthru, reverse, passthru),
        "B": (passthru, passthru, reverse),
        "R": (reverse, passthru, passthru),
        # N, P, C all use identity functions (extra logic applied later)
        "N": (passthru, passthru, passthru),
        "P": (passthru, passthru, passthru),
        "C": (passthru, passthru, passthru),
    }

    if rev_mode not in funcs:
        raise ValueError(f"Unsupported rev mode: {rev_mode}")

    return funcs[rev_mode]


def _quality_ok(seq_qual: str, asize: int, minqual: float) -> bool:
    """
    Legacy-style quality check: convert ASCII to uint8, subtract 35,
    then require the mean of each anchor region to be >= minqual.
    """
    if len(seq_qual) < asize * 2:
        return False

    # fromstring is deprecated; frombuffer reproduces same behaviour.
    arr = np.frombuffer(seq_qual.encode("ascii"), dtype=np.uint8) - 35
    left_ok = float(arr[:asize].mean()) >= minqual
    right_ok = float(arr[-asize:].mean()) >= minqual
    return left_ok and right_ok


def _emit_anchor_pair(
    read_name: str,
    seq: str,
    qual: str,
    asize: int,
    read_f,
    A_f,
    B_f,
    rev_mode: str,
):
    """
    Emit the A/B anchors for a single read to stdout in FASTQ format.
    This mirrors the final print-block in the legacy script.
    """
    # Apply read-level transform (A/B/R/N/P/C)
    seq, qual = read_f(seq), read_f(qual)

    # Optional permutation ("P" control)
    # The full permutation pool & burn-in is handled outside this function;
    # here we just assume seq/qual are already permuted if needed.

    # Optional full reverse-complement ("C")
    if rev_mode == "C":
        seq = rev_comp(seq)
        qual = qual[::-1]

    # Emit A anchor
    sys.stdout.write(f"@{read_name}_A__{seq}\n")
    sys.stdout.write(A_f(seq[:asize]) + "\n")
    sys.stdout.write("+\n")
    sys.stdout.write(A_f(qual[:asize]) + "\n")

    # Emit B anchor
    sys.stdout.write(f"@{read_name}_B\n")
    sys.stdout.write(B_f(seq[-asize:]) + "\n")
    sys.stdout.write("+\n")
    sys.stdout.write(B_f(qual[-asize:]) + "\n")


def _drive_reads_from_bam(
    path: Path,
) -> Iterator[pysam.AlignedSegment]:
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            yield read


def _drive_reads_from_reads_file(
    path: Path,
) -> Iterator[Tuple[str, str]]:
    # sites.reads: "name<TAB>sequence"
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            name, seq = line.split("\t")
            yield name.replace(" ", "_"), seq


def _drive_reads_from_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    with path.open() as fh:
        for name, seq in fasta_chunks(fh):
            yield name.replace(" ", "_"), seq


def _drive_reads_from_fastq(path: Path) -> Iterator[Tuple[str, str, str]]:
    with path.open() as fh:
        for q in qfa_chunks(fh):
            yield q.name.replace(" ", "_"), q.seq, q.qual


@click.command()
@click.argument("input_path", type=click.Path(exists=True, readable=True))
@click.option(
    "-a",
    "--anchor",
    "asize",
    type=int,
    default=20,
    show_default=True,
    help="Anchor size.",
)
@click.option(
    "-q",
    "--min-qual",
    "minqual",
    type=float,
    default=5.0,
    show_default=True,
    help="Min avg. quality on both anchors.",
)
@click.option(
    "-r",
    "--rev",
    "rev_mode",
    type=click.Choice(["A", "B", "R", "C", "N", "P"]),
    default="N",
    show_default=True,
    help="Reverse / permutation mode (as in legacy unmapped2anchors.py).",
)
@click.option(
    "-R",
    "--reads",
    is_flag=True,
    help="Input is sites.reads (name<TAB>seq) instead of BAM.",
)
@click.option(
    "-F",
    "--fasta",
    is_flag=True,
    help="Input is FASTA instead of BAM.",
)
@click.option(
    "-Q",
    "--fastq",
    is_flag=True,
    help="Input is FASTQ instead of BAM.",
)
def main(
    input_path: str,
    asize: int,
    minqual: float,
    rev_mode: str,
    reads: bool,
    fasta: bool,
    fastq: bool,
):
    """
    Split unmapped reads into A/B anchors (Python 3 port of unmapped2anchors.py).
    """
    path = Path(input_path)

    # Ensure mutually exclusive input modes
    mode_flags = [reads, fasta, fastq]
    if sum(bool(x) for x in mode_flags) > 1:
        raise click.UsageError("Use at most one of --reads/--fasta/--fastq")

    # Build transform functions
    read_f, A_f, B_f = _build_funcs(rev_mode)

    # Permutation control state (for rev_mode == "P")
    perm_A = []
    perm_I = []
    perm_B = []
    perm_burn_in = []
    N_perm = 100

    def handle_seq(name: str, seq: str, qual: str):
        nonlocal perm_A, perm_I, perm_B, perm_burn_in

        if len(seq) < asize * 2:
            return

        # quality filter (legacy behaviour)
        if not _quality_ok(qual, asize=asize, minqual=minqual):
            return

        # Permutation control: collect burn-in, then randomise A/I/B independently
        if rev_mode == "P":
            perm_A.append((seq[:asize], qual[:asize]))
            perm_B.append((seq[-asize:], qual[-asize:]))
            perm_I.append((seq[asize:-asize], qual[asize:-asize]))

            if len(perm_burn_in) < N_perm:
                perm_burn_in.append((name, seq, qual))
                return

            import random

            A_seq, A_qual = perm_A.pop(random.randrange(len(perm_A)))
            B_seq, B_qual = perm_B.pop(random.randrange(len(perm_B)))
            I_seq, I_qual = perm_I.pop(random.randrange(len(perm_I)))

            seq = A_seq + I_seq + B_seq
            qual = A_qual + I_qual + B_qual

        _emit_anchor_pair(
            read_name=name,
            seq=seq,
            qual=qual,
            asize=asize,
            read_f=read_f,
            A_f=A_f,
            B_f=B_f,
            rev_mode=rev_mode,
        )

    # --- Drive reads depending on mode ------------------------------------ #

    if reads:
        for idx, (name, seq) in enumerate(_drive_reads_from_reads_file(path), start=1):
            qual = "b" * len(seq)
            handle_seq(f"{name}_{idx}", seq, qual)
    elif fasta:
        for idx, (name, seq) in enumerate(_drive_reads_from_fasta(path), start=1):
            qual = "b" * len(seq)
            handle_seq(f"{name}_{idx}", seq, qual)
    elif fastq:
        for idx, (name, seq, qual) in enumerate(_drive_reads_from_fastq(path), start=1):
            handle_seq(f"{name}_{idx}", seq, qual)
    else:
        # BAM mode: take only unmapped reads (as in legacy)
        for read in _drive_reads_from_bam(path):
            if read.is_unmapped:
                seq = read.query_sequence or ""
                qual_arr = read.query_qualities
                if seq is None or qual_arr is None:
                    continue
                # Convert numeric qualities back to ASCII string; assume Phred+33
                qual = "".join(chr(q + 33) for q in qual_arr)
                handle_seq(read.query_name, seq, qual)

    # After main pass, emit perm_burn_in if in P-mode
    for name, seq, qual in perm_burn_in:
        handle_seq(name, seq, qual)


if __name__ == "__main__":
    main()
