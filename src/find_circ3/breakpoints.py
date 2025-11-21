from __future__ import annotations

from dataclasses import dataclass
from typing import List, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .engine import AnchorHit, FindCircConfig, GenomeAccessor


@dataclass
class BreakpointHit:
    """
    Representation of a single candidate breakpoint between two anchors.

    Modernised counterpart of the tuples returned by the legacy
    find_breakpoints() in find_circ.py. The legacy tuples had:

        (dist, ov, strandmatch, rnd, chrom, start, end, signal, sense)

    Here we keep deterministic fields and let the caller decide what to do
    with ties.
    """

    dist: int           # edit distance between internal read and spliced flanks
    overlap: int        # anchor overlap in base pairs
    strandmatch: str    # "MATCH", "MISMATCH" or "NA"
    chrom: str
    start: int          # genomic start coordinate of the junction
    end: int            # genomic end coordinate of the junction
    signal: str         # 4bp splice signal, e.g. "GTAG", "GCAG"
    sense: str          # '+' or '-' (annotated splice sense)


def mismatches(a: str, b: str) -> int:
    """
    Count the number of mismatching characters between two sequences.

    Sequences are truncated to the same length before comparison.
    """
    if not a or not b:
        return max(len(a), len(b))

    n = min(len(a), len(b))
    a_arr = np.frombuffer(a[:n].encode("ascii"), dtype="uint8")
    b_arr = np.frombuffer(b[:n].encode("ascii"), dtype="uint8")
    return int((a_arr != b_arr).sum())


def _aligned_length(cigar: str) -> int:
    """
    Return the number of reference bases consumed by an alignment CIGAR.

    Counts M, D, N, =, X as consuming reference positions. This mirrors
    pysam's reference_length behaviour and the legacy use of A.aend.
    """
    if not cigar:
        return 0

    length = 0
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        if not num:
            continue
        n = int(num)
        if ch in ("M", "D", "N", "=", "X"):
            length += n
        num = ""
    return length


def _rev_comp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def find_breakpoints(
    A: "AnchorHit",
    B: "AnchorHit",
    read_seq: str,
    chrom: str,
    cfg: "FindCircConfig",
    genome: "GenomeAccessor",
) -> List[BreakpointHit]:
    """
    Estimate splice breakpoints between two anchors A,B supported by a
    full read sequence.

    This is a faithful but slightly simplified port of the legacy
    find_breakpoints() in find_circ.py:

      - reconstruct an internal read segment between both anchors,
      - take genomic flanks around A and B from the reference genome,
      - for each possible breakpoint position x inside the internal segment:
          * compute edit distance between the "spliced" genome sequence and
            the internal read,
          * estimate anchor overlap using the same margin/l formula as
            the original,
          * derive a 4bp splice signal and a coarse strandmatch label,
      - return all breakpoints that are tied by (dist, overlap), and if
        requested, also by strandmatch.
    """
    asize = cfg.anchor_size
    # Defaults analogous to legacy options.margin and options.maxdist
    margin = getattr(cfg, "margin", max(1, asize // 4))
    maxdist = getattr(cfg, "max_mismatches", 2)
    strandpref = getattr(cfg, "strandpref", False)

    L = len(read_seq)
    if L == 0:
        return []

    # eff_a = options.asize - margin in the legacy code
    eff_a = asize - margin
    if eff_a <= 0 or 2 * eff_a >= L:
        # Not enough internal sequence to do a meaningful search
        return []

    internal = read_seq[eff_a : L - eff_a].upper()
    l = len(internal)
    if l <= 0:
        return []

    flank_len = l + 2  # as in the original script

    # Compute anchor end position for A based on CIGAR.
    a_ref_len = _aligned_length(A.cigar)
    if a_ref_len <= 0:
        return []

    # In our engine, A.pos is 1-based; the last aligned base is:
    aend = A.pos + a_ref_len - 1

    # Genomic windows around A and B – analogue of:
    #  A_flank = genome.get(chrom, A.aend-margin, A.aend-margin + flank, '+')
    #  B_flank = genome.get(chrom, B.pos - flank+margin, B.pos+margin, '+')
    # but using get_seq(chrom, start, end) with 1-based inclusive coords.
    a_start = max(1, aend - margin)
    a_end = a_start + flank_len - 1

    b_start = B.pos - flank_len + margin
    if b_start < 1:
        b_start = 1
    b_end = b_start + flank_len - 1

    try:
        A_flank = genome.get_seq(chrom, a_start, a_end).upper()
        B_flank = genome.get_seq(chrom, b_start, b_end).upper()
    except Exception:
        return []

    # Ensure we have enough context; if not, shrink to the available length.
    max_l = min(len(A_flank), len(B_flank)) - 2
    if max_l <= 0:
        return []
    if l > max_l:
        l = max_l
        internal = internal[:l]
        flank_len = l + 2
        A_flank = A_flank[:flank_len]
        B_flank = B_flank[:flank_len]

    hits: List[BreakpointHit] = []
    canonical = {"GTAG", "GCAG", "ATAC"}

    for x in range(l + 1):
        # Splice genome flanks at position x:
        #   first x bases from A_flank,
        #   then everything from B_flank after the donor site (x+2).
        spliced = A_flank[:x] + B_flank[x + 2 :]
        dist = mismatches(spliced, internal)

        if dist > maxdist:
            continue

        # Anchor overlap as in the original code.
        ov = 0
        if x < margin:
            ov = margin - x
        if l - x < margin:
            ov = margin - (l - x)

        gt = A_flank[x : x + 2]
        ag = B_flank[x : x + 2]
        gtag = (gt + ag).upper()
        rc_gtag = _rev_comp(gtag)

        # Determine signal and sense based on canonical patterns.
        if gtag in canonical:
            signal = gtag
            sense = "+"
            strandmatch = "MATCH"
        elif rc_gtag in canonical:
            signal = gtag
            sense = "-"
            strandmatch = "MATCH"
        else:
            signal = gtag
            sense = "+"
            strandmatch = "MISMATCH"

        # For now, keep these equal to the anchor geometry; the caller
        # uses its own junction coords for BED.
        start = min(A.pos, B.pos)
        end = max(A.pos, B.pos)

        hits.append(
            BreakpointHit(
                dist=dist,
                overlap=ov,
                strandmatch=strandmatch,
                chrom=chrom,
                start=start,
                end=end,
                signal=signal,
                sense=sense,
            )
        )

    if not hits:
        return []

    # Hits are sorted, with low edit distance beating low anchor overlap.
    # This matches "hits = sorted(hits)" in the legacy implementation.
    hits.sort(key=lambda h: (h.dist, h.overlap, h.strandmatch))

    best = hits[0]

    if strandpref:
        # Exploit strand information to break ties if requested:
        # keep only those tied in dist, overlap *and* strandmatch.
        ties = [
            h
            for h in hits
            if (h.dist == best.dist)
            and (h.overlap == best.overlap)
            and (h.strandmatch == best.strandmatch)
        ]
    else:
        # Default: tie only on edit distance and overlap.
        ties = [
            h
            for h in hits
            if (h.dist == best.dist) and (h.overlap == best.overlap)
        ]

    return ties
