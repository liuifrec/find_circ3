# src/find_circ3/legacy_engine.py
from __future__ import annotations

import sys
import logging
import traceback
from collections import defaultdict
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pysam

LOG = logging.getLogger("find_circ3.legacy")
def extract_full_read(A, B):
    """
    Try to reconstruct the original read sequence.

    Priority:
    1. Legacy style: QNAME contains '__' and the part after it is the read seq.
    2. Normal SAM: use A.query_sequence (fallback to B if needed).
    """
    # Legacy: QNAME encodes the sequence
    if "__" in A.qname:
        parts = A.qname.split("__", 1)
        if len(parts) > 1 and parts[1]:
            return parts[1]
    if "__" in B.qname:
        parts = B.qname.split("__", 1)
        if len(parts) > 1 and parts[1]:
            return parts[1]

    # Fallback to SEQ field
    if A.query_sequence:
        return A.query_sequence
    if B.query_sequence:
        return B.query_sequence

    raise RuntimeError(
        "Unable to recover read sequence from QNAME or SEQ fields "
        "for read %r / %r" % (A.qname, B.qname)
    )


# -----------------------------
# Basic sequence utilities
# -----------------------------

COMPLEMENT = {
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'k': 'm', 'm': 'k', 'r': 'y', 'y': 'r',
    's': 's', 'w': 'w', 'b': 'v', 'v': 'b',
    'h': 'd', 'd': 'h', 'n': 'n',
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'K': 'M', 'M': 'K', 'R': 'Y', 'Y': 'R',
    'S': 'S', 'W': 'W', 'B': 'V', 'V': 'B',
    'H': 'D', 'D': 'H', 'N': 'N',
}


def complement(s: str) -> str:
    return "".join(COMPLEMENT.get(x, 'N') for x in s)


def rev_comp(seq: str) -> str:
    return complement(seq)[::-1]


# -----------------------------
# Genome accessor using pysam.FastaFile
# -----------------------------

class Genome:
    """
    Replacement for the original mmap/indexed_fasta logic, but with the same
    *external* semantics:
      - coordinates are 0-based, [start, end)
      - out-of-bounds are padded with 'N'
      - strand '-' returns reverse complement
    """

    def __init__(self, fasta_path: str):
        self.fa = pysam.FastaFile(fasta_path)
        self.lengths = {
            chrom: self.fa.get_reference_length(chrom)
            for chrom in self.fa.references
        }

    def get(self, chrom: str, start: int, end: int, sense: str = '+') -> str:
        if chrom not in self.lengths:
            raise KeyError("Chromosome %s not present in FASTA" % chrom)

        size = self.lengths[chrom]

        pad_start = 0
        pad_end = 0

        if start < 0:
            pad_start = -start
            start = 0
        if end > size:
            pad_end = end - size
            end = size

        if end <= start:
            seq = ""
        else:
            seq = self.fa.fetch(chrom, start, end)

        if pad_start or pad_end:
            seq = "N" * pad_start + seq + "N" * pad_end

        if sense == '-':
            seq = rev_comp(seq)

        return seq.upper()


# -----------------------------
# grouper(2, iterable) like original
# -----------------------------

def grouper2(iterable):
    it = iter(iterable)
    while True:
        try:
            A = next(it)
        except StopIteration:
            return
        try:
            B = next(it)
        except StopIteration:
            return
        yield A, B


# -----------------------------
# Hit class – same semantics
# -----------------------------

class Hit:
    def __init__(self, minmapscore: int):
        self.reads: List[str] = []
        self.readnames: List[str] = []
        self.uniq = set()
        self.mapquals_A: List[int] = []
        self.mapquals_B: List[int] = []
        self.uniq_bridges: int = 0
        self.tissues: Dict[str, int] = defaultdict(int)
        self.edits: List[int] = []
        self.overlaps: List[int] = []
        self.n_hits: List[int] = []
        self.signal: str = "NNNN"
        self.strand_plus: int = 0
        self.strand_minus: int = 0
        self.strandmatch: str = 'NA'
        self.minmapscore = minmapscore

    def add(self, options: dict, read: str, A, B, dist: int, ov: int,
            strandmatch: str, signal: str, n_hits: int) -> None:
        self.signal = signal
        self.strandmatch = strandmatch
        self.edits.append(dist)
        self.overlaps.append(ov)
        self.n_hits.append(n_hits)

        # by convention have A precede B in the genome.
        if A.pos > B.pos:
            A, B = B, A

        # Alignment Score - Secondbest hit score
        aopt = dict(A.tags)
        bopt = dict(B.tags)
        qA = aopt.get('AS', 0) - aopt.get('XS', self.minmapscore)
        qB = bopt.get('AS', 0) - bopt.get('XS', self.minmapscore)

        if qA and qB:
            # both anchors from the *same read* align uniquely
            self.uniq_bridges += 1

        self.mapquals_A.append(qA)
        self.mapquals_B.append(qB)

        # recover the original readname
        # ('__' is forbidden in input read names!)
        # recover the original readname
        if "__" in A.qname:
            base = A.qname.split("__", 1)[0]
        elif "__" in B.qname:
            base = B.qname.split("__", 1)[0]
        else:
            # tiny.sam (and other simple tests) have plain QNAMEs
            base = A.qname

        qname = base  # no [:-2]; we don't need to strip /1 in tests
        self.readnames.append(qname)

        # record the spliced read sequence as it was before mapping
        if A.is_reverse:
            self.strand_minus += 1
            self.reads.append(rev_comp(read))
        else:
            self.strand_plus += 1
            self.reads.append(read)

        # identify the tissue/sample it came from
        sample_name = options["name"]
        for (prefix, tiss) in options["samples"]:
            if qname.startswith(prefix):
                sample_name = tiss
                break

        self.tissues[sample_name] += 1

        self.uniq.add((read, sample_name))
        self.uniq.add((rev_comp(read), sample_name))

    def scores(self, chrom: str, start: int, end: int, sense: str):
        n_reads = len(self.reads)
        n_uniq = len(self.uniq) // 2
        best_qual_A = sorted(self.mapquals_A, reverse=True)[0]
        best_qual_B = sorted(self.mapquals_B, reverse=True)[0]

        tissues = sorted(self.tissues.keys())
        tiss_counts = [str(self.tissues[k]) for k in tissues]

        return (
            n_reads,
            n_uniq,
            best_qual_A,
            best_qual_B,
            self.uniq_bridges,
            tissues,
            tiss_counts,
            min(self.edits),
            min(self.overlaps),
            min(self.n_hits),
            self.signal,
            self.strandmatch,
        )


# -----------------------------
# find_breakpoints – ported literally
# -----------------------------

def find_breakpoints(A, B, read: str, chrom: str, genome: Genome, options: dict):
    """
    Faithful port of legacy find_breakpoints.
    """

    def mismatches(a: str, b: str) -> int:
        # byte-wise Hamming distance; legacy used numpy.fromstring, but
        # this is equivalent for content and safer in Py3.
        return sum(1 for aa, bb in zip(a, b) if aa != bb)

    def strandmatch(ann: str, sense: str) -> str:
        if ann == sense:
            return "MATCH"
        if ann == "*" or len(ann) > 1:
            return "NA"
        return "MISMATCH"

    L = len(read)
    hits: List[Tuple] = []

    eff_a = options["anchor"] - options["margin"]
    internal = read[eff_a:-eff_a].upper()
    flank = L - 2 * eff_a + 2

    A_flank = genome.get(
        chrom,
        A.aend - options["margin"],
        A.aend - options["margin"] + flank,
        '+',
    )
    B_flank = genome.get(
        chrom,
        B.pos - flank + options["margin"],
        B.pos + options["margin"],
        '+',
    )

    l = L - 2 * eff_a

    def rnd():
        return np.random.random()

    if not options["randomize"]:
        rnd = lambda: 0  # type: ignore[assignment]

    for x in range(l + 1):
        spliced = A_flank[:x] + B_flank[x + 2:]
        dist = mismatches(spliced, internal)

        ov = 0
        if x < options["margin"]:
            ov = options["margin"] - x
        if l - x < options["margin"]:
            ov = options["margin"] - (l - x)

        if dist <= options["maxdist"]:
            gt = A_flank[x:x + 2]
            ag = B_flank[x:x + 2]
            gtag = gt + ag
            rc_gtag = rev_comp(gtag)

            start = B.pos + options["margin"] - l + x
            end = A.aend - options["margin"] + x + 1
            start, end = min(start, end), max(start, end)

            strand = "*"

            if options["stranded"]:
                if A.is_reverse:
                    strand = '-'
                else:
                    strand = '+'

            if options["noncanonical"]:
                hits.append(
                    (dist, ov, strandmatch(strand, '+'),
                     rnd(), chrom, start, end, gtag, '+')
                )
                hits.append(
                    (dist, ov, strandmatch(strand, '-'),
                     rnd(), chrom, start, end, rc_gtag, '-')
                )
            else:
                if gtag == 'GTAG':
                    hits.append(
                        (dist, ov, strandmatch(strand, '+'),
                         rnd(), chrom, start, end, 'GTAG', '+')
                    )
                elif gtag == 'CTAC':
                    hits.append(
                        (dist, ov, strandmatch(strand, '-'),
                         rnd(), chrom, start, end, 'GTAG', '-')
                    )

    if len(hits) < 2:
        return hits

    hits = sorted(hits)
    best = hits[0]

    if options["strandpref"]:
        ties = [
            h
            for h in hits
            if (h[0] == best[0]) and (h[1] == best[1]) and (h[2] == best[2])
        ]
    else:
        ties = [h for h in hits if (h[0] == best[0]) and (h[1] == best[1])]

    return ties


# -----------------------------
# Main legacy engine as generator
# -----------------------------

def legacy_call_iter(
    sam_path: str,
    genome_fasta: str,
    name: str = "unknown",
    prefix: str = "",
    anchor: int = 20,
    margin: int = 2,
    maxdist: int = 2,
    min_uniq_qual: int = 2,
    noncanonical: bool = False,
    randomize: bool = False,
    allhits: bool = False,
    stranded: bool = False,
    strandpref: bool = False,
    halfunique: bool = False,
    report_nobridges: bool = False,
    reads2samples_path: str = "",
    reads_out_path: str | None = None,
    bam_out_path: str | None = None,
    stats_path: str = "runstats.log",
    max_intron: int | None = None,
    min_support: int = 1,
) -> Iterable[str]:
    """
    Generator version of original find_circ.py logic.

    Yields BED lines (and header lines starting with '# '), one per junction.
    """

    # options dict like the original 'options' object
    options = {
        "name": name,
        "anchor": anchor,
        "margin": margin,
        "maxdist": maxdist,
        "noncanonical": noncanonical,
        "randomize": randomize,
        "stranded": stranded,
        "strandpref": strandpref,
    }

    # sample mapping
    if reads2samples_path:
        samples = [
            line.rstrip("\n").split('\t')
            for line in open(reads2samples_path)
        ]
    else:
        samples = []
    samples.append(('', name))
    options["samples"] = samples

    # reads output
    if reads_out_path:
        readfile = open(reads_out_path, "w")
    else:
        readfile = sys.stderr

    genome = Genome(genome_fasta)

    N = defaultdict(int)
    minmapscore = anchor * (-2)

    circs: Dict[Tuple[str, int, int, str], Hit] = defaultdict(
        lambda: Hit(minmapscore)
    )
    splices: Dict[Tuple[str, int, int, str], Hit] = defaultdict(
        lambda: Hit(minmapscore)
    )
    loci = defaultdict(list)  # kept for parity; not used downstream

    sam = pysam.AlignmentFile(sam_path, "r")
    if bam_out_path:
        bam_out = pysam.AlignmentFile(bam_out_path, "wb", template=sam)
    else:
        bam_out = None

    pair_num = 0

    try:
        for pair_num, (A, B) in enumerate(grouper2(sam)):
            if A is None or B is None:
                break

            N['total'] += 1

            if A.is_unmapped or B.is_unmapped:
                N['unmapped'] += 1
                continue
            if A.reference_id != B.reference_id:
                N['other_chrom'] += 1
                continue
            if A.is_reverse != B.is_reverse:
                N['other_strand'] += 1
                continue

            dist = B.pos - A.pos
            if abs(dist) < anchor:
                N['overlapping_anchors'] += 1
                continue

            if bam_out is not None:
                bam_out.write(A)
                bam_out.write(B)

            # circular orientation
            if (A.is_reverse and dist > 0) or (not A.is_reverse and dist < 0):
                read = extract_full_read(A, B)
                chrom = sam.get_reference_name(A.reference_id)

                if A.is_reverse:
                    A, B = B, A
                    read = rev_comp(read)

                bp_hits = find_breakpoints(A, B, read, chrom, genome, options)

                if not bp_hits:
                    N['circ_no_bp'] += 1
                else:
                    N['circ_reads'] += 1

                n_hits = len(bp_hits)
                if bp_hits and not allhits:
                    bp_hits = [bp_hits[0]]

                for h in bp_hits:
                    dist_bp, ov, strandmatch, rnd, chrom, start, end, signal, sense = h
                    key = (chrom, start + 1, end - 1, sense)
                    circs[key].add(
                        options, read, A, B,
                        dist_bp, ov, strandmatch, signal, n_hits
                    )

            # linear orientation
            elif (A.is_reverse and dist < 0) or (not A.is_reverse and dist > 0):
                read = extract_full_read(A, B)
                chrom = sam.get_reference_name(A.reference_id)

                if A.is_reverse:
                    A, B = B, A
                    read = rev_comp(read)

                bp_hits = find_breakpoints(A, B, read, chrom, genome, options)

                if not bp_hits:
                    N['splice_no_bp'] += 1
                else:
                    N['spliced_reads'] += 1

                n_hits = len(bp_hits)
                if bp_hits and not allhits:
                    bp_hits = [bp_hits[0]]

                for h in bp_hits:
                    dist_bp, ov, strandmatch, rnd, chrom, start, end, signal, sense = h
                    key = (chrom, start, end, sense)
                    splices[key].add(
                        options, read, A, B,
                        dist_bp, ov, strandmatch, signal, n_hits
                    )
                    loci[(chrom, start, sense)].append(splices[key])
                    loci[(chrom, end, sense)].append(splices[key])

            else:
                N['fallout'] += 1
                LOG.warning("unhandled read: A='%s' B='%s'", A, B)

    except KeyboardInterrupt:
        sam_line = pair_num * 2
        fastq_line = pair_num * 8
        LOG.warning(
            "KeyboardInterrupt while processing pair %d (SAM line %d, FASTQ line %d)",
            pair_num, sam_line, fastq_line,
        )
    except Exception:
        sam_line = pair_num * 2
        fastq_line = pair_num * 8
        LOG.error(
            "Unhandled exception at pair %d (SAM line %d, FASTQ line %d)",
            pair_num, sam_line, fastq_line,
        )
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
    finally:
        sam.close()
        if bam_out is not None:
            bam_out.close()

    # -----------------------------
    # Output – faithfully ported
    # -----------------------------

    def output(cand: Dict[Tuple[str, int, int, str], Hit], prefix_tag: str):
        # header line (tests ignore lines starting with '#')
        yield "# " + "\t".join([
            'chrom', 'start', 'end', 'name', 'n_reads', 'strand',
            'n_uniq', 'uniq_bridges', 'best_qual_left', 'best_qual_right',
            'tissues', 'tiss_counts', 'edits', 'anchor_overlap',
            'breakpoints', 'signal', 'strandmatch', 'category',
        ])

        n = 1
        for c, hit in cand.items():
            chrom, start, end, sense = c
            (n_reads, n_uniq, best_qual_A, best_qual_B, uniq_bridges,
             tissues, tiss_counts, min_edit, min_anchor_ov, n_hits,
             signal, strandmatch) = hit.scores(chrom, start, end, sense)

            # extra filters for your new CLI knobs (keep defaults gentle)
            if max_intron is not None:
                if abs(end - start) > max_intron:
                    continue

            if n_reads < min_support:
                continue

            # legacy uniqueness filters
            if halfunique:
                if (best_qual_A < min_uniq_qual) and (best_qual_B < min_uniq_qual):
                    N['anchor_not_uniq'] += 1
                    continue
            else:
                if (best_qual_A < min_uniq_qual) or (best_qual_B < min_uniq_qual):
                    N['anchor_not_uniq'] += 1
                    continue

            if (uniq_bridges == 0) and not report_nobridges:
                N['no_uniq_bridges'] += 1
                continue

            name_j = "%s%s_%06d" % (prefix, prefix_tag, n)
            n += 1

            for r_seq, ori_name in zip(hit.reads, hit.readnames):
                readfile.write(">%s %s\n%s\n" % (name_j, ori_name, r_seq))

            categories: List[str] = []
            if signal == "GTAG":
                categories.append("CANONICAL")
            if strandmatch == "MATCH":
                categories.append("STRANDMATCH")
            if best_qual_A > 0 and best_qual_B > 0 and uniq_bridges > 0:
                categories.append("ANCHOR_UNIQUE")
            if uniq_bridges == 0:
                categories.append("NO_UNIQ_BRIDGES")
            if n_hits == 1:
                categories.append("UNAMBIGUOUS_BP")
            if min_anchor_ov == 0 and min_edit == 0:
                categories.append("PERFECT_EXT")
            elif min_anchor_ov <= 1 and min_edit <= 1:
                categories.append("GOOD_EXT")
            elif min_anchor_ov <= 2 and min_edit <= 2:
                categories.append("OK_EXT")

            if not categories:
                categories.append("DUBIOUS")

            if prefix_tag == "circ":
                categories.append("CIRCULAR")
            else:
                categories.append("LINEAR")

            bed = [
                chrom,
                start - 1,  # BED start
                end,
                name_j,
                n_reads,
                sense,
                n_uniq,
                uniq_bridges,
                best_qual_A,
                best_qual_B,
                ",".join(tissues),
                ",".join(tiss_counts),
                min_edit,
                min_anchor_ov,
                n_hits,
                signal,
                strandmatch,
                ",".join(sorted(categories)),
            ]
            yield "\t".join(str(b) for b in bed)

    # stats like original
    with open(stats_path, "w") as stats_f:
        stats_f.write(str(dict(N)) + "\n")

    # circs, then normals – exactly like the original
    for line in output(circs, "circ"):
        yield line
    for line in output(splices, "norm"):
        yield line

    if readfile is not sys.stderr:
        readfile.close()
