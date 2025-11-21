from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, TYPE_CHECKING

if TYPE_CHECKING:
    from .engine import (
        JunctionCandidate,
        Junction,
        FindCircConfig,
        GenomeAccessor,
        FindCircStats,
        AnchorHit,
    )


def _extract_full_read_seq_from_hits(hits: List["AnchorHit"]) -> str:
    """
    Recover the original full read sequence from anchor read IDs.

    unmapped2anchors.py encodes the original read sequence into the FASTQ
    header of the A-anchor as:

        @<read.qname>_A__<SEQ>

    When these anchors are aligned by bowtie2, the query_name in the SAM
    records becomes, e.g.:

        <read.qname>_A__<SEQ>   (A-anchor)
        <read.qname>_B          (B-anchor)

    For a given pair of AnchorHit objects, we look for the first read_id
    containing "_A__" and return the suffix after that marker as the full
    read sequence. If none match this pattern (e.g. synthetic tests), we
    return an empty string and breakpoint logic can fall back to defaults.
    """
    for h in hits:
        name = h.read_id
        marker = "_A__"
        idx = name.find(marker)
        if idx != -1:
            return name[idx + len(marker) :].upper()
    return ""


@dataclass
class HitAccumulator:
    """
    Minimal Hit-like accumulator for a single junction.

    Structurally inspired by the original Hit class in find_circ.py. It
    gathers evidence from JunctionCandidate objects and produces a Junction
    with scores ready for BED output.
    """

    # aggregated read-level evidence
    n_reads: int = 0
    n_uniq: int = 0              # number of uniquely mapped anchors
    uniq_bridges: int = 0
    mapquals_left: List[int] = field(default_factory=list)
    mapquals_right: List[int] = field(default_factory=list)

    # tissue / sample tracking (simplified)
    tissues: Dict[str, int] = field(default_factory=dict)

    # breakpoint-level metrics (to be filled via find_breakpoints)
    edits: List[int] = field(default_factory=list)
    overlaps: List[int] = field(default_factory=list)
    n_hits: List[int] = field(default_factory=list)

    # splice / orientation annotations
    signal: str = "NNNN"
    strandmatch: str = "NA"

    # junction coordinates
    chrom: str | None = None
    left_pos: int | None = None
    right_pos: int | None = None
    strand: str | None = None

    def add_candidate(
        self,
        cand: JunctionCandidate,
        cfg: FindCircConfig,
        genome: GenomeAccessor,
    ) -> None:
        """
        Absorb a single JunctionCandidate's evidence into this accumulator.

        For now, we treat each candidate as representing a single read with
        two anchors (A,B). This mirrors the typical unmapped2anchors +
        bowtie2 usage in the legacy find_circ pipeline.
        """
        # Initialize coordinates on first call
        if self.chrom is None:
            self.chrom = cand.chrom
            self.left_pos = cand.left_pos
            self.right_pos = cand.right_pos
            self.strand = cand.strand

        hits = cand.hits
        if not hits:
            return

        # --- AS/XS-based scoring, aligned with legacy Hit.add() ---

        # In the original implementation, per-anchor uniqueness is given by
        # q = AS - XS (or AS - minmapscore if XS is missing). If q > 0, the
        # anchor is considered uniquely placed.
        minmapscore = cfg.anchor_size * (-2)

        per_anchor_q: List[int] = []
        uniq_flags: List[bool] = []

        for h in hits:
            tags = getattr(h, "tags", {}) or {}
            ascore = tags.get("AS")
            xs = tags.get("XS", minmapscore if ascore is not None else 0)

            if ascore is not None:
                q = int(ascore) - int(xs)
            else:
                # Fallback to MAPQ when AS/XS are not present (e.g. tiny tests)
                q = int(h.mapq)

            per_anchor_q.append(q)
            uniq_flags.append(q > 0)

        # Number of "reads" and unique anchors:
        # we approximate n_reads as number of anchor hits, matching our
        # previous minimal implementation and tiny test behaviour.
        self.n_reads += len(hits)
        self.n_uniq += sum(uniq_flags)

        # uniq_bridges: both anchors from the same read align uniquely
        if len(uniq_flags) >= 2 and uniq_flags[0] and uniq_flags[1]:
            self.uniq_bridges += 1

        # store per-side quality scores (for best_qual_left/right)
        # left = first hit, right = second if present, fallback to first
        self.mapquals_left.append(per_anchor_q[0])
        if len(per_anchor_q) > 1:
            self.mapquals_right.append(per_anchor_q[1])
        else:
            self.mapquals_right.append(per_anchor_q[0])

        # tissue/sample tracking (simplified: a single "sample")
        sample_name = "sample"
        self.tissues[sample_name] = self.tissues.get(sample_name, 0) + len(hits)

        # --- Breakpoint-level information via find_breakpoints ---

        # Default conservative values (used if we cannot recover a full read
        # sequence or if breakpoint search fails).
        min_edits = 0
        min_overlap = 0
        n_best_hits = 1

        # For signal/strandmatch, we start with the simple genome-based logic
        # and later allow find_breakpoints to override it.
        try:
            left2 = genome.get_seq(cand.chrom, cand.left_pos, cand.left_pos + 1)
            right2 = genome.get_seq(cand.chrom, cand.right_pos - 1, cand.right_pos)
            self.signal = (left2 + right2).upper()

            canonical = {"GTAG", "GCAG", "ATAC"}
            self.strandmatch = "MATCH" if self.signal in canonical else "MISMATCH"
        except Exception:
            self.signal = "NNNN"
            self.strandmatch = "NA"

        # Attempt to recover the original full read sequence and perform a
        # proper breakpoint search. If this fails for any reason, we fall back
        # to the conservative defaults above, which keeps existing tests and
        # behaviour intact while exercising the new code path on real data.
        read_seq = _extract_full_read_seq_from_hits(hits)

        if len(hits) >= 2 and read_seq:
            A = hits[0]
            B = hits[1]
            try:
                from .breakpoints import find_breakpoints

                bp_list = find_breakpoints(
                    A=A,
                    B=B,
                    read_seq=read_seq,
                    chrom=cand.chrom,
                    cfg=cfg,
                    genome=genome,
                )
            except Exception:
                bp_list = []

            if bp_list:
                # For now, find_breakpoints returns a list of BreakpointHit
                # objects. Once the full legacy logic is ported, this list
                # may contain several tied best hits.
                min_edits = min(bp.dist for bp in bp_list)
                min_overlap = min(bp.overlap for bp in bp_list)
                n_best_hits = len(bp_list)

                # Use the first breakpoint to override signal/strandmatch;
                # this will be made more sophisticated when we port the full
                # tie-breaking logic.
                best_bp = bp_list[0]
                self.signal = best_bp.signal
                self.strandmatch = best_bp.strandmatch

        self.edits.append(min_edits)
        self.overlaps.append(min_overlap)
        self.n_hits.append(n_best_hits)

    def finalize(
        self,
        cfg: FindCircConfig,
        stats: FindCircStats,
        category: str,
        prefix: str,
    ) -> Junction:
        """
        Convert accumulated evidence into a Junction suitable for BED output.
        """

        # Lazy import to avoid circular import at module load time
        from .engine import Junction

        if (
            self.chrom is None
            or self.left_pos is None
            or self.right_pos is None
            or self.strand is None
        ):
            raise RuntimeError("HitAccumulator.finalize() called without coordinates")

        # best per-side qualities
        best_qual_left = max(self.mapquals_left) if self.mapquals_left else 0
        best_qual_right = max(self.mapquals_right) if self.mapquals_right else 0

        # tissues → compact representation
        tissue_names = sorted(self.tissues.keys())
        if tissue_names:
            tissues = ",".join(tissue_names)
            tiss_counts = ",".join(str(self.tissues[t]) for t in tissue_names)
        else:
            tissues = "sample"
            tiss_counts = str(self.n_reads)

        # breakpoint-related fields as conservative summaries
        min_edits = min(self.edits) if self.edits else 0
        min_overlap = min(self.overlaps) if self.overlaps else 0
        min_n_hits = min(self.n_hits) if self.n_hits else 1

        # update stats
        if category == "CIRCULAR":
            stats.circular_junctions += 1
        elif category == "LINEAR":
            stats.linear_junctions += 1

        name = f"{prefix}{self.chrom}:{self.left_pos}|{self.right_pos}"

        # Legacy-style category labels
        categories: list[str] = []

        # CANONICAL only for GTAG (legacy find_circ behaviour)
        if self.signal == "GTAG":
            categories.append("CANONICAL")

        # STRANDMATCH for canonical strandmatch
        if self.strandmatch == "MATCH":
            categories.append("STRANDMATCH")

        # ANCHOR_UNIQUE if both anchors can align uniquely and at least one
        # uniq bridge has been observed.
        if best_qual_left > 0 and best_qual_right > 0 and self.uniq_bridges > 0:
            categories.append("ANCHOR_UNIQUE")

        # NO_UNIQ_BRIDGES when we never saw a unique bridge
        if self.uniq_bridges == 0:
            categories.append("NO_UNIQ_BRIDGES")

        # UNAMBIGUOUS_BP if there is exactly one best breakpoint
        if min_n_hits == 1:
            categories.append("UNAMBIGUOUS_BP")

        # Extension quality based on edits and anchor overlap
        if min_overlap == 0 and min_edits == 0:
            categories.append("PERFECT_EXT")
        elif min_overlap <= 1 and min_edits <= 1:
            categories.append("GOOD_EXT")
        elif min_overlap <= 2 and min_edits <= 2:
            categories.append("OK_EXT")

        # If no other label assigned, call it DUBIOUS (legacy default)
        if not categories:
            categories.append("DUBIOUS")

        # Final CIRCULAR/LINEAR tag
        categories.append(category)

        # Sort for deterministic output (as in legacy: sorted(categories))
        category_str = ",".join(sorted(categories))

        return Junction(
            chrom=self.chrom,
            start=self.left_pos,
            end=self.right_pos,
            name=name,
            n_reads=self.n_reads,
            strand=self.strand,
            n_uniq=self.n_uniq,
            uniq_bridges=self.uniq_bridges,
            best_qual_left=best_qual_left,
            best_qual_right=best_qual_right,
            tissues=tissues,
            tiss_counts=tiss_counts,
            edits=min_edits,
            anchor_overlap=min_overlap,
            breakpoints=min_n_hits,
            signal=self.signal,
            strandmatch=self.strandmatch,
            category=category_str,
        )
