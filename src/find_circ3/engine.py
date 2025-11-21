from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional, List, Literal

import logging
import pysam

from .hit_accumulator import HitAccumulator

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public types and containers
# ---------------------------------------------------------------------------

JunctionType = Literal["circular", "linear", "other"]


@dataclass
class AnchorHit:
    """
    Representation of a single anchor alignment (one hit from bowtie2).

    This is roughly analogous to one SAM alignment record for an anchor read
    in the legacy find_circ implementation.
    """

    read_id: str
    chrom: str
    pos: int          # 1-based leftmost position on the reference (SAM-style)
    strand: str       # "+" or "-"
    cigar: str
    mapq: int
    is_spliced: bool  # True if the CIGAR contains an 'N' (intron)

    # New: keep full tag dictionary and reverse flag (for future use)
    tags: dict = field(default_factory=dict)
    is_reverse: bool = False


@dataclass
class JunctionCandidate:
    """
    A potential splice junction supported by one read (or read pair),
    derived from one or more AnchorHit objects.

    This is an intermediate representation between raw alignments and the
    final Junction objects that are written to BED.
    """

    read_id: str
    chrom: str
    left_pos: int
    right_pos: int
    strand: str

    hits: List[AnchorHit]

    # early classification before full scoring
    tentative_type: JunctionType = "other"


@dataclass
class FindCircConfig:
    """
    Configuration for a find_circ3 run.

    This is the high-level, "Pythonic" interface that mirrors the options
    of the original find_circ.py script, but without exposing CLI concerns.
    """

    anchors_fastq: Path
    genome: Path
    sample_name: str = "unknown"
    prefix: str = ""
    min_uniq_qual: int = 2
    anchor_size: int = 20
    stats_path: Optional[Path] = None
    reads_path: Optional[Path] = None


@dataclass
class FindCircStats:
    """
    Run-time statistics, roughly mirroring the counters that find_circ
    writes into its log file (circ_no_bp, lin_no_bp, etc.).
    """

    total_reads: int = 0
    processed_anchors: int = 0
    candidate_junctions: int = 0
    circular_junctions: int = 0
    linear_junctions: int = 0

    # generic container for any extra counters we might port later
    extra: dict[str, float] = field(default_factory=dict)


@dataclass
class Junction:
    """
    Representation of a single splice junction as produced by find_circ.

    This mirrors the 18-column BED-like output format described in the
    original find_circ README (columns 1–18).
    """

    chrom: str
    start: int
    end: int
    name: str
    n_reads: int
    strand: str

    n_uniq: int
    uniq_bridges: int
    best_qual_left: int
    best_qual_right: int

    tissues: str
    tiss_counts: str

    edits: int
    anchor_overlap: int
    breakpoints: int
    signal: str
    strandmatch: str
    category: str


# ---------------------------------------------------------------------------
# Genome accessor
# ---------------------------------------------------------------------------


class GenomeAccessor:
    """
    Minimal genome accessor abstraction.

    In the original find_circ implementation, this wraps an mmap'ed FASTA
    and provides fast substring access. Here we start with a simple pysam-
    based accessor for a single FASTA file.

    Notes
    -----
    - We currently support only a single multi-FASTA file.
    - Folder-of-chromosomes support can be added later if needed.
    """

    def __init__(self, genome_path: Path):
        self.genome_path = genome_path
        if genome_path.is_file():
            # Single FASTA file
            self._fasta = pysam.FastaFile(str(genome_path))
        else:
            # Placeholder for folder-based references
            self._fasta = None
            raise NotImplementedError(
                "Folder-based genome references are not implemented yet."
            )

    def get_seq(self, chrom: str, start: int, end: int) -> str:
        """
        Return genomic sequence for [start, end] (1-based, inclusive).

        pysam's fetch uses 0-based, end-exclusive coordinates, so we convert
        appropriately.
        """
        if self._fasta is None:
            raise NotImplementedError("GenomeAccessor has no FASTA loaded.")
        # Convert 1-based inclusive → 0-based [start-1, end)
        return self._fasta.fetch(chrom, start - 1, end)


def _open_genome(genome: Path) -> GenomeAccessor:
    return GenomeAccessor(genome)


# ---------------------------------------------------------------------------
# Anchor reading and candidate detection
# ---------------------------------------------------------------------------


def _iter_anchors(anchors_path: Path) -> Iterator[AnchorHit]:
    """
    Iterate over anchor alignments for the find_circ3 engine.

    This opens a SAM/BAM (bowtie2 or bwa output) using pysam and yields
    one AnchorHit per alignment.
    """

    suffix = anchors_path.suffix.lower()
    if suffix == ".sam":
        mode = "r"   # text SAM
    else:
        mode = "rb"  # BAM/CRAM-style

    with pysam.AlignmentFile(str(anchors_path), mode) as af:
        for aln in af.fetch(until_eof=True):
            if aln.is_unmapped:
                continue

            read_id = aln.query_name
            chrom = af.get_reference_name(aln.reference_id)
            pos = int(aln.reference_start) + 1  # 0-based → 1-based
            strand = "-" if aln.is_reverse else "+"
            cigar = aln.cigarstring or ""
            mapq = int(aln.mapping_quality)
            is_spliced = "N" in cigar

            tags = dict(aln.tags)

            yield AnchorHit(
                read_id=read_id,
                chrom=chrom,
                pos=pos,
                strand=strand,
                cigar=cigar,
                mapq=mapq,
                is_spliced=is_spliced,
                tags=tags,
                is_reverse=aln.is_reverse,
            )



def _candidate_from_pair(
    A: AnchorHit,
    B: AnchorHit,
    cfg: FindCircConfig,
) -> Optional[JunctionCandidate]:
    """
    Construct a JunctionCandidate from a single (A,B) anchor pair.

    This mirrors the core geometry logic from the original find_circ.py:
      - anchors must be on the same chromosome and strand,
      - pairs closer than anchor_size are discarded as overlapping anchors,
      - reversed orientation pairs → circular candidates,
      - sequential orientation pairs → linear candidates.
    """

    # Require same chromosome
    if A.chrom != B.chrom:
        return None

    # Require same strand
    if A.strand != B.strand:
        return None

    # Distance in genomic coordinates
    dist = B.pos - A.pos
    if abs(dist) < cfg.anchor_size:
        # overlapping anchors; do not form a candidate
        return None

    A_is_reverse = (A.strand == "-")

    # Legacy logic:
    #   if (A.is_reverse and dist > 0) or (not A.is_reverse and dist < 0):
    #       → circular
    #   elif (A.is_reverse and dist < 0) or (not A.is_reverse and dist > 0):
    #       → linear
    if (A_is_reverse and dist > 0) or (not A_is_reverse and dist < 0):
        tentative_type: JunctionType = "circular"
    elif (A_is_reverse and dist < 0) or (not A_is_reverse and dist > 0):
        tentative_type = "linear"
    else:
        return None

    left_pos = min(A.pos, B.pos)
    right_pos = max(A.pos, B.pos)

    return JunctionCandidate(
        read_id=A.read_id,
        chrom=A.chrom,
        left_pos=left_pos,
        right_pos=right_pos,
        strand=A.strand,
        hits=[A, B],
        tentative_type=tentative_type,
    )


def _detect_junction_candidates(
    anchor_stream: Iterator[AnchorHit],
    genome_accessor: GenomeAccessor,
    cfg: FindCircConfig,
    stats: FindCircStats,
) -> Iterator[JunctionCandidate]:
    """
    Core junction candidate discovery step.

    This mirrors the original find_circ behaviour more closely: it consumes
    the anchor_stream in pairs (A,B), assuming that each pair of consecutive
    alignments belongs together (as produced by unmapped2anchors + bowtie2).
    """

    buffer: AnchorHit | None = None

    for hit in anchor_stream:
        if buffer is None:
            # start a new pair
            buffer = hit
            continue

        # We have a full pair (A,B)
        A = buffer
        B = hit
        buffer = None

        # Stats: one "read" worth of anchors processed
        stats.total_reads += 1
        stats.processed_anchors += 2

        cand = _candidate_from_pair(A, B, cfg)
        if cand is not None:
            stats.candidate_junctions += 1
            yield cand

    # If there's an odd number of alignments, we silently drop the last one.


# ---------------------------------------------------------------------------
# Scoring & filtering via HitAccumulator
# ---------------------------------------------------------------------------


def _score_and_filter_junctions(
    candidate_junctions: Iterator[JunctionCandidate],
    genome_accessor: GenomeAccessor,
    cfg: FindCircConfig,
    stats: FindCircStats,
) -> Iterator[Junction]:
    """
    Scoring and filtering of candidate junctions.

    Delegates most of the aggregation to HitAccumulator, which is
    structurally aligned with the original Hit class from find_circ.py.
    """

    type_to_category = {
        "circular": "CIRCULAR",
        "linear": "LINEAR",
        "other": "OTHER",
    }

    for cand in candidate_junctions:
        if cand.tentative_type not in ("circular", "linear"):
            continue

        category = type_to_category.get(cand.tentative_type, "OTHER")
        if category == "OTHER":
            continue

        acc = HitAccumulator()
        acc.add_candidate(cand, cfg, genome_accessor)
        j = acc.finalize(cfg, stats, category=category, prefix=cfg.prefix)
        yield j


# ---------------------------------------------------------------------------
# Stats / supporting reads / BED output
# ---------------------------------------------------------------------------


def _write_stats_if_requested(stats: FindCircStats, cfg: FindCircConfig) -> None:
    """
    Write numeric run statistics if cfg.stats_path is set.

    The exact contents will be modelled after the original find_circ log.
    For now this is a placeholder.
    """
    if cfg.stats_path is None:
        return
    logger.debug("Stats writing is not implemented yet; skipping.")


def _write_supporting_reads_if_requested(cfg: FindCircConfig) -> None:
    """
    Write supporting reads to a separate file if requested.

    In the original script, supporting reads are written to stderr by
    default, or to a file if an option is provided.
    """
    if cfg.reads_path is None:
        return
    logger.debug("Supporting reads writing is not implemented yet; skipping.")


def _junctions_to_bed(
    final_junctions: Iterator[Junction],
    cfg: FindCircConfig,
) -> Iterator[str]:
    """
    Convert final junction objects to BED-like lines.

    We mirror the column structure of the original find_circ
    splice_sites.bed output (18 columns).
    """

    for j in final_junctions:
        fields = [
            j.chrom,
            str(j.start),
            str(j.end),
            j.name,
            str(j.n_reads),
            j.strand,
            str(j.n_uniq),
            str(j.uniq_bridges),
            str(j.best_qual_left),
            str(j.best_qual_right),
            j.tissues,
            j.tiss_counts,
            str(j.edits),
            str(j.anchor_overlap),
            str(j.breakpoints),
            j.signal,
            j.strandmatch,
            j.category,
        ]
        yield "\t".join(fields)


# ---------------------------------------------------------------------------
# Public engine API
# ---------------------------------------------------------------------------


def run_find_circ(
    anchors_fastq: Path,
    genome: Path,
    sample_name: str = "unknown",
    prefix: str = "",
    min_uniq_qual: int = 2,
    anchor_size: int = 20,
    stats_path: Optional[Path] = None,
    reads_path: Optional[Path] = None,
) -> Iterable[str]:
    """
    Core find_circ3 engine.

    Designed to reproduce the behaviour of the original find_circ script,
    but with a modern, testable API.
    """
    cfg = FindCircConfig(
        anchors_fastq=anchors_fastq,
        genome=genome,
        sample_name=sample_name,
        prefix=prefix,
        min_uniq_qual=min_uniq_qual,
        anchor_size=anchor_size,
        stats_path=stats_path,
        reads_path=reads_path,
    )
    stats = FindCircStats()

    logger.info("find_circ3 starting: %s", cfg)

    genome_accessor = _open_genome(cfg.genome)
    anchor_stream = _iter_anchors(cfg.anchors_fastq)
    candidate_junctions = _detect_junction_candidates(
        anchor_stream,
        genome_accessor,
        cfg,
        stats,
    )
    final_junctions = _score_and_filter_junctions(
        candidate_junctions,
        genome_accessor,
        cfg,
        stats,
    )

    _write_stats_if_requested(stats, cfg)
    _write_supporting_reads_if_requested(cfg)

    for bed_line in _junctions_to_bed(final_junctions, cfg):
        yield bed_line
