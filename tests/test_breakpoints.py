from pathlib import Path

from find_circ3.breakpoints import find_breakpoints
from find_circ3.engine import AnchorHit, FindCircConfig


class FakeGenome:
    def __init__(self, seq: str):
        self.seq = seq.upper()

    def get_seq(self, chrom: str, start: int, end: int) -> str:
        # 1-based inclusive coordinates
        start = max(1, start)
        end = min(len(self.seq), end)
        if end < start:
            return ""
        return self.seq[start - 1 : end]


def test_breakpoints_edits_and_overlap_simple():
    # Construct a toy genome of length 40 and plant specific 4bp flanks
    # around the A and B anchors such that the breakpoint search finds
    # exactly one best hit with dist=0 and overlap=0.
    seq = list("N" * 40)

    # A_flank will be genome[12:15] (1-based, len 4): "AGTT"
    a_flank = "AGTT"
    for i, ch in enumerate(a_flank, start=12 - 1):
        seq[i] = ch

    # B_flank will be genome[27:30] (1-based, len 4): "CCAG"
    b_flank = "CAGG"

    for i, ch in enumerate(b_flank, start=27 - 1):
        seq[i] = ch

    genome_seq = "".join(seq)
    genome = FakeGenome(genome_seq)

    # Read sequence: anchor_size=4, margin defaults to 1 (anchor_size//4),
    # so eff_a = 3 and internal = read_seq[3:len-3] = "AG" for this read.
    read_seq = "AAAAGGAA"  # length 8 -> internal "AG"
    asize = 4

    A = AnchorHit(
        read_id="read_A__" + read_seq,
        chrom="chr1",
        pos=10,         # arbitrary, but consistent with flank placement
        strand="+",
        cigar="4M",
        mapq=60,
        is_spliced=False,
        tags={},
        is_reverse=False,
    )

    B = AnchorHit(
        read_id="read_B",
        chrom="chr1",
        pos=30,
        strand="+",
        cigar="4M",
        mapq=60,
        is_spliced=False,
        tags={},
        is_reverse=False,
    )

    cfg = FindCircConfig(
        anchors_fastq=Path("dummy.fastq"),
        genome=Path("dummy.fa"),
        sample_name="test",
        prefix="",
        min_uniq_qual=2,
        anchor_size=asize,
        stats_path=None,
        reads_path=None,
    )

    bps = find_breakpoints(
        A=A,
        B=B,
        read_seq=read_seq,
        chrom="chr1",
        cfg=cfg,
        genome=genome,
    )

    assert bps, "Expected at least one breakpoint candidate"

    # Our construction ensures a unique best candidate with:
    #  - dist == 0 (perfect spliced vs internal)
    #  - overlap == 0 (x in the central region of the window)
    bp = bps[0]
    assert bp.dist == 0
    assert bp.overlap == 0

    # At the chosen x, flanks produce a canonical GT-AG signal.
    assert bp.signal == "GTAG"
    assert bp.strandmatch == "MATCH"
