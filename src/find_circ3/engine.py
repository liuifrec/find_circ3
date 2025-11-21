from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional, TextIO
from collections import defaultdict
from typing import Dict
from typing import Optional

import logging

from typing import List, Literal
import pysam

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public configuration & stats containers
# ---------------------------------------------------------------------------
JunctionType = Literal["circular", "linear", "other"]
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

    # Decide circular vs linear based on the same conditions as the legacy code:
    #
    #   if (A.is_reverse and dist > 0) or (not A.is_reverse and dist < 0):
    #       # reversed orientation → circRNA
    #   elif (A.is_reverse and dist < 0) or (not A.is_reverse and dist > 0):
    #       # sequential → linear
    #
    if (A_is_reverse and dist > 0) or (not A_is_reverse and dist < 0):
        tentative_type: JunctionType = "circular"
    elif (A_is_reverse and dist < 0) or (not A_is_reverse and dist > 0):
        tentative_type = "linear"
    else:
        # fallout case; we skip it
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



    
@dataclass
class AnchorHit:
    """
    Representation of a single anchor alignment (one hit from bowtie2).

    In the legacy implementation this corresponds roughly to one SAM
    alignment record for an anchor read.
    """

    read_id: str
    chrom: str
    pos: int          # 1-based leftmost position on the reference (SAM-style)
    strand: str       # "+" or "-"
    cigar: str
    mapq: int
    is_spliced: bool  # True if the CIGAR contains an 'N' (intron)


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
    Run-time statistics, roughly mirroring the counters that find_circ2
    writes into its log file (circ_no_bp, lin_no_bp, etc.).

    For now this is just a placeholder; we'll fill individual fields while
    porting the logic from the legacy implementation.
    """

    total_reads: int = 0
    processed_anchors: int = 0
    candidate_junctions: int = 0
    circular_junctions: int = 0
    linear_junctions: int = 0

    # We keep a generic dictionary for extra counters we may want to track
    # one-to-one with the original implementation.
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
# Public API
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

    This function is intended to reproduce the behaviour of the original
    `find_circ.py` script from find_circ2, but with a modern, testable API.

    Parameters
    ----------
    anchors_fastq:
        FASTQ file containing the anchor reads (output of an
        unmapped2anchors-like step).
    genome:
        Reference genome FASTA file or folder with per-chromosome FASTAs.
    sample_name:
        Name of the sample/tissue. Used for naming junctions and stats.
    prefix:
        Optional prefix to prepend to each junction identifier.
    min_uniq_qual:
        Minimal uniqueness score for anchor alignments.
    anchor_size:
        Anchor size used when generating anchors.
    stats_path:
        Optional path to write numeric run statistics.
    reads_path:
        Optional path to write supporting reads instead of stderr.

    Returns
    -------
    Iterable[str]
        An iterable over BED-like lines (strings), each representing a
        linear or circular junction. This should match the semantics of
        the original find_circ2 output (e.g. CIRCULAR / LINEAR tags).
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

    # NOTE:
    # The following high-level steps are deliberately structured to map to
    # sections of the original find_circ.py implementation. When we port
    # the algorithm, each TODO below will become a concrete, well-scoped
    # function.

    # 1) Open genome accessor (ported from GenomeAccessor in find_circ2)
    genome_accessor = _open_genome(cfg.genome)

    # 2) Stream anchor alignments / reads
    #    In the original code, this effectively reads a SAM stream from
    #    bowtie2. In find_circ3 we will *eventually* either:
    #      - read SAM/BAM from stdin, or
    #      - operate on a pre-generated "anchors SAM" file.
    #
    #    For now, this is kept as a placeholder iterator.
    anchor_stream = _iter_anchors(cfg.anchors_fastq)

    # 3) Generate candidate junctions from anchor pairs
    #    This corresponds to the "breaking anchor" logic and the core
    #    back-splice detection code in find_circ2.
    candidate_junctions = _detect_junction_candidates(
        anchor_stream,
        genome_accessor,
        cfg,
        stats,
    )

    # 4) Score, filter, and classify junctions (CIRCULAR vs LINEAR)
    final_junctions = _score_and_filter_junctions(
    candidate_junctions,
    genome_accessor,
    cfg,
    stats,
    )


    # 5) Optionally write stats and supporting reads
    _write_stats_if_requested(stats, cfg)
    _write_supporting_reads_if_requested(cfg)

    # 6) Convert final junctions to BED-like strings and yield them.
    for bed_line in _junctions_to_bed(final_junctions, cfg):
        yield bed_line


# ---------------------------------------------------------------------------
# Internal stubs to be filled during porting
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
    # TODO: in the original code, this builds a GenomeAccessor with mmap.
    return GenomeAccessor(genome)



def _iter_anchors(anchors_path: Path) -> Iterator[AnchorHit]:
    """
    Iterate over anchor alignments for the find_circ3 engine.

    This opens a SAM/BAM (bowtie2 or bwa output) using pysam and yields
    one AnchorHit per alignment.

    For future regression tests, this will typically read something like
    'legacy_anchors.sam' derived from the original find_circ test_data.
    """

    # Auto-detect SAM vs BAM based on suffix. This keeps things simple for now.
    suffix = anchors_path.suffix.lower()
    if suffix == ".sam":
        mode = "r"   # text SAM
    else:
        # default to BAM/CRAM-style binary input
        mode = "rb"

    with pysam.AlignmentFile(str(anchors_path), mode) as af:
        for aln in af.fetch(until_eof=True):
            # skip unmapped reads; they don't contribute anchor hits
            if aln.is_unmapped:
                continue

            read_id = aln.query_name
            chrom = af.get_reference_name(aln.reference_id)
            # SAM is 0-based start → convert to 1-based
            pos = int(aln.reference_start) + 1
            strand = "-" if aln.is_reverse else "+"
            cigar = aln.cigarstring or ""
            mapq = int(aln.mapping_quality)
            is_spliced = "N" in cigar

            yield AnchorHit(
                read_id=read_id,
                chrom=chrom,
                pos=pos,
                strand=strand,
                cigar=cigar,
                mapq=mapq,
                is_spliced=is_spliced,
            )




def _detect_junction_candidates(
    anchor_stream: Iterator[AnchorHit],
    genome_accessor: GenomeAccessor,
    cfg: FindCircConfig,
    stats: FindCircStats,
) -> Iterator[JunctionCandidate]:
    """
    Core junction candidate discovery step.

    This implementation mirrors the original find_circ behavior more closely:
    it consumes the anchor_stream in pairs (A,B), assuming that each pair of
    consecutive alignments belongs together (as produced by unmapped2anchors
    + bowtie2).

    For each (A,B) pair, we apply the same-chromosome, same-strand, distance
    and orientation rules to decide whether to emit a circular or linear
    JunctionCandidate.
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
    # This matches the behavior of grouper(2, sam) in the original code.


def _score_and_filter_junctions(
    candidate_junctions: Iterator[JunctionCandidate],
    genome_accessor: GenomeAccessor,
    cfg: FindCircConfig,
    stats: FindCircStats,
) -> Iterator[Junction]:
    """
    Scoring and filtering of candidate junctions.

    This is a minimal, first-pass implementation inspired by the legacy
    Hit.add / Hit.scores logic from find_circ:

      - aggregate per-candidate read support,
      - compute simple uniqueness / bridge metrics from MAPQ,
      - derive a 4bp splice signal from the genome,
      - assign a category ("CIRCULAR"/"LINEAR") based on tentative_type.

    It does *not* yet replicate full breakpoint refinement or all of the
    original filters; those can be added incrementally.
    """

    # Simple helper to map tentative_type → final BED category
    type_to_category = {
        "circular": "CIRCULAR",
        "linear": "LINEAR",
        "other": "OTHER",
    }

    for cand in candidate_junctions:
        # For now, we only emit circular/linear candidates
        if cand.tentative_type not in ("circular", "linear"):
            continue

        hits = cand.hits
        if not hits:
            continue

        # Core support metrics (Hit-like)
        n_reads = len(hits)
        # count "unique" hits based on mapq threshold
        uniq_flags = [h.mapq >= cfg.min_uniq_qual for h in hits]
        n_uniq = sum(uniq_flags)
        # uniq_bridges: both anchors from same read align uniquely
        uniq_bridges = 1 if all(uniq_flags) and len(hits) >= 2 else 0

        # Best qualities for A/B (take first two hits if available)
        best_qual_left = hits[0].mapq
        best_qual_right = hits[1].mapq if len(hits) > 1 else hits[0].mapq

        # Tissues / counts: we don't track tissues yet; use a simple placeholder
        tissues = "sample"
        tiss_counts = str(n_reads)

        # Edits / overlaps / n_hits: not implemented yet; placeholders
        edits = 0
        anchor_overlap = 0
        breakpoints = 0

        # Splice signal and strandmatch using genome
        signal = "NNNN"
        strandmatch = "NA"

        try:
            # Very simple 2+2 bp signal: 2bp around left and right positions.
            # In a more faithful port, we'd use the refined breakpoints;
            # here we approximate with the candidate's left/right coords.
            left2 = genome_accessor.get_seq(cand.chrom, cand.left_pos, cand.left_pos + 1)
            right2 = genome_accessor.get_seq(cand.chrom, cand.right_pos - 1, cand.right_pos)
            signal = (left2 + right2).upper()

            # Canonical splice motifs (GT-AG, GC-AG, AT-AC); this is a simplification.
            canonical = {"GTAG", "GCAG", "ATAC"}
            strandmatch = "MATCH" if signal in canonical else "MISMATCH"
        except Exception:
            # If genome access fails, keep defaults
            signal = "NNNN"
            strandmatch = "NA"

        category = type_to_category.get(cand.tentative_type, "OTHER")

        # Update high-level stats
        if category == "CIRCULAR":
            stats.circular_junctions += 1
        elif category == "LINEAR":
            stats.linear_junctions += 1

        # Construct a name similar to find_circ IDs
        name = f"{cfg.prefix}{cand.chrom}:{cand.left_pos}|{cand.right_pos}"

        j = Junction(
            chrom=cand.chrom,
            start=cand.left_pos,
            end=cand.right_pos,
            name=name,
            n_reads=n_reads,
            strand=cand.strand,
            n_uniq=n_uniq,
            uniq_bridges=uniq_bridges,
            best_qual_left=best_qual_left,
            best_qual_right=best_qual_right,
            tissues=tissues,
            tiss_counts=tiss_counts,
            edits=edits,
            anchor_overlap=anchor_overlap,
            breakpoints=breakpoints,
            signal=signal,
            strandmatch=strandmatch,
            category=category,
        )
        yield j


def _write_stats_if_requested(stats: FindCircStats, cfg: FindCircConfig) -> None:
    """
    Write numeric run statistics if cfg.stats_path is set.

    The exact contents will be modelled after the original find_circ2 log
    (circ_no_bp, lin_no_bp, etc.). For now, this is a placeholder.
    """
    if cfg.stats_path is None:
        return
    # TODO: implement once stats fields are properly populated.
    logger.debug("Stats writing is not implemented yet; skipping.")


def _write_supporting_reads_if_requested(cfg: FindCircConfig) -> None:
    """
    Write supporting reads to a separate file if requested.

    In the original script, supporting reads are written to stderr by
    default, or to a file if an option is provided.
    """
    if cfg.reads_path is None:
        return
    # TODO: implement once we have the representation of supporting reads.
    logger.debug("Supporting reads writing is not implemented yet; skipping.")


def _junctions_to_bed(
    final_junctions: Iterator[Junction],
    cfg: FindCircConfig,
) -> Iterator[str]:
    """
    Convert final junction objects to BED-like lines.

    We mirror the column structure and semantics of the original
    find_circ splice_sites.bed output (18 columns).
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
