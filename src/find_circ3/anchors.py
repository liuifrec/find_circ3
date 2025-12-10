from __future__ import annotations

from typing import Iterable, TextIO
import statistics

import click
import pysam

# Standard Sanger / Illumina 1.8+ Phred offset
QUAL_OFFSET = 33

# Legacy unmapped2anchors.py used: nq = numpy.fromstring(qual, dtype=uint8) - 35
LEGACY_OFFSET = 35


def _legacy_scores(qual: str) -> list[int]:
    """
    Convert a FASTQ quality string into 'legacy' pseudo-Phred scores,
    matching:

        nq = numpy.fromstring(qual, dtype=uint8) - 35

    in the original Py2 script.
    """
    return [ord(ch) - LEGACY_OFFSET for ch in qual]



def _phred_to_qual_string(qualities: Iterable[int]) -> str:
    """
    Convert an iterable of Phred scores (ints, numpy arrays, or array.array)
    to a FASTQ quality string.
    """
    return "".join(chr(int(q) + QUAL_OFFSET) for q in qualities)

def _emit_anchor_fastq(
    out: TextIO,
    header: str,            # full header *without* leading '@'
    seq: str,
    quals: Iterable[int],
) -> None:
    qual_str = _phred_to_qual_string(quals)
    out.write(f"@{header}\n")
    out.write(seq + "\n+\n" + qual_str + "\n")


def _anchors_from_bam(path: str, anchor_size: int, min_qual: int, out: TextIO) -> int:
    n = 0
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            if not read.is_unmapped:
                continue

            seq = read.query_sequence
            quals = read.query_qualities
            qual_str = read.qual  # raw ASCII qualities, needed for legacy-style filter

            if seq is None or quals is None or qual_str is None:
                continue
            if len(seq) < 2 * anchor_size or len(qual_str) < 2 * anchor_size:
                continue

            # --- Legacy-style quality filter (mean over each anchor) ---
            legacy = _legacy_scores(qual_str)
            left_mean = statistics.mean(legacy[:anchor_size])
            right_mean = statistics.mean(legacy[-anchor_size:])

            if left_mean < min_qual or right_mean < min_qual:
                # "read is junk" in unmapped2anchors.py
                continue

            # Anchor sequences and Phred scores (for writing)
            left_seq  = seq[:anchor_size]
            right_seq = seq[-anchor_size:]

            left_qual  = quals[:anchor_size]
            right_qual = quals[-anchor_size:]

            qname = read.query_name
            full = seq  # full read sequence, for the A-header

            # Py2 behaviour:
            #   A-header: <qname>_A__<FULL_READ_SEQ>
            #   B-header: <qname>_B
            _emit_anchor_fastq(
                out=out,
                header=f"{qname}_A__{full}",
                seq=left_seq,
                quals=left_qual,
            )
            _emit_anchor_fastq(
                out=out,
                header=f"{qname}_B",
                seq=right_seq,
                quals=right_qual,
            )
            n += 2
    return n

def _anchors_from_fastq(path: str, anchor_size: int, min_qual: int, out: TextIO) -> int:
    n = 0
    with open(path, "r", encoding="ascii") as fh:
        idx = 0
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not (seq and plus and qual):
                break

            header = header.rstrip("\n\r")
            seq = seq.rstrip("\n\r")
            qual = qual.rstrip("\n\r")

            if not header.startswith("@"):
                continue

            # Py2 qfa_chunks + handle_read:
            #   read.qname = "%s_%d" % (name, N)
            idx += 1
            name = header[1:].replace(" ", "_")
            base_name = f"{name}_{idx}"

            if len(seq) < 2 * anchor_size or len(qual) < 2 * anchor_size:
                continue

            # --- Legacy-style quality filter (mean over each anchor) ---
            legacy = _legacy_scores(qual)
            left_mean = statistics.mean(legacy[:anchor_size])
            right_mean = statistics.mean(legacy[-anchor_size:])

            if left_mean < min_qual or right_mean < min_qual:
                continue

            # For writing, we still use Phred ints with QUAL_OFFSET=33
            phred = [ord(ch) - QUAL_OFFSET for ch in qual]

            left_seq  = seq[:anchor_size]
            right_seq = seq[-anchor_size:]
            left_qual  = phred[:anchor_size]
            right_qual = phred[-anchor_size:]

            full = seq
            _emit_anchor_fastq(
                out=out,
                header=f"{base_name}_A__{full}",
                seq=left_seq,
                quals=left_qual,
            )
            _emit_anchor_fastq(
                out=out,
                header=f"{base_name}_B",
                seq=right_seq,
                quals=right_qual,
            )
            n += 2
    return n


@click.command()
@click.argument(
    "input_path",
    type=click.Path(exists=True, dir_okay=False, readable=True),
)
@click.option(
    "--fastq",
    is_flag=True,
    default=False,
    help="Interpret INPUT_PATH as FASTQ instead of BAM.",
)
@click.option(
    "--anchor",
    "anchor_size",
    type=int,
    required=True,
    help="Anchor size in bases.",
)
@click.option(
    "--min-qual",
    type=int,
    default=35,
    show_default=True,
    help="Minimum 'legacy' quality threshold (see docs).",
)
def main(input_path: str, fastq: bool, anchor_size: int, min_qual: int) -> None:
    """
    Emit anchor reads in FASTQ format from a BAM (default) or FASTQ (--fastq).

    For each input read that passes filtering, two anchors are written:
    the leftmost and rightmost `--anchor` bases.
    """
    import sys

    if fastq:
        _anchors_from_fastq(input_path, anchor_size, min_qual, sys.stdout)
    else:
        _anchors_from_bam(input_path, anchor_size, min_qual, sys.stdout)
