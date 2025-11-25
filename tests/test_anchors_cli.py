# tests/test_anchors_cli.py

from pathlib import Path

import pysam
from click.testing import CliRunner

from find_circ3.anchors import main as anchors_main


def _make_unmapped_bam(tmp_path: Path) -> Path:
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 50}]}
    bam_path = tmp_path / "unmapped.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out_bam:
        aln = pysam.AlignedSegment()
        aln.query_name = "read1"
        aln.query_sequence = "ACGTACGTACGTACGT"
        aln.flag = 4  # unmapped
        aln.query_qualities = pysam.qualitystring_to_array("BBBBBBBBBBBBBBBB")
        out_bam.write(aln)
    return bam_path


def test_anchors_from_bam(tmp_path):
    bam = _make_unmapped_bam(tmp_path)
    runner = CliRunner()
    result = runner.invoke(
        anchors_main,
        [str(bam), "--anchor", "4", "--min-qual", "0"],
    )
    assert result.exit_code == 0, result.output

    lines = [ln.strip() for ln in result.output.splitlines() if ln.strip()]
    # We expect 8 lines: 4 for A, 4 for B
    assert len(lines) == 8

    # Header names
    assert lines[0].startswith("@read1_A__")
    assert lines[4].startswith("@read1_B")


def test_anchors_from_fastq(tmp_path):
    fq = tmp_path / "reads.fastq"
    fq.write_text(
        "@r1\n"
        "ACGTACGT\n"
        "+\n"
        "BBBBBBBB\n"
    )

    runner = CliRunner()
    result = runner.invoke(
        anchors_main,
        [str(fq), "--fastq", "--anchor", "4", "--min-qual", "0"],
    )
    assert result.exit_code == 0, result.output

    lines = [ln.strip() for ln in result.output.splitlines() if ln.strip()]
    assert len(lines) == 8
    assert lines[0].startswith("@r1_1_A__")
    assert lines[4].startswith("@r1_1_B")
