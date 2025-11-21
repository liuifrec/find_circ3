from pathlib import Path
import gzip

from click.testing import CliRunner

from find_circ3.cli import main as cli_main


def test_cli_tiny_sam_matches_expected():
    data_dir = Path("tests/test_against_legacy/data")

    tiny_sam = data_dir / "tiny.sam"
    tiny_genome = data_dir / "tiny_genome.fa"
    expected_bed = data_dir / "tiny_expected.bed"

    assert tiny_sam.exists(), "tiny.sam missing"
    assert tiny_genome.exists(), "tiny_genome.fa missing"
    assert expected_bed.exists(), "tiny_expected.bed missing"

    runner = CliRunner()
    result = runner.invoke(
    cli_main,
        [
            "call",
            str(tiny_sam),
            "--genome",
            str(tiny_genome),
            "--name",
            "test",
            "--prefix",
            "test_",
            "--anchor",
            "5",           # <-- add this
        ],
        )

    # Engine should run without NotImplementedError now
    assert result.exit_code == 0, result.output

    # Collect non-empty, non-comment lines from CLI output
    out_lines = [
        line.strip()
        for line in result.output.splitlines()
        if line.strip() and not line.startswith("#")
    ]

    with expected_bed.open() as f:
        exp_lines = [
            line.strip()
            for line in f
            if line.strip() and not line.startswith("#")
        ]

    assert out_lines == exp_lines


def test_legacy_files_exist_and_are_nonempty():
    data_dir = Path("tests/test_against_legacy/data")

    anchors_gz = data_dir / "legacy_anchors.fastq.gz"
    anchors_fq = data_dir / "legacy_anchors.fastq"
    bed = data_dir / "legacy_splice_sites.bed"

    # Accept either gzipped or plain FASTQ
    if anchors_gz.exists():
        anchors = anchors_gz
        opener = lambda p: gzip.open(p, "rt")
    else:
        anchors = anchors_fq
        opener = lambda p: open(p, "r")

    assert anchors.exists(), "Legacy anchors FASTQ is missing"
    assert bed.exists(), "Legacy splice sites BED is missing"

    # anchors should have some reads (FASTQ line count > 0)
    with opener(anchors) as f:
        n_lines = sum(1 for _ in f)
    assert n_lines > 0, "Legacy anchors FASTQ is empty"

    # BED should have at least one non-comment line
    with bed.open() as f:
        entries = [line for line in f if line.strip() and not line.startswith("#")]
    assert len(entries) > 0, "Legacy splice_sites BED has no entries"
def test_cli_cdr1as_emits_circular_junctions():
    data_dir = Path("tests/test_against_legacy/data")

    cdr1as_sam = data_dir / "cdr1as_anchors.sam"
    cdr1as_genome = data_dir / "CDR1as_locus.fa"

    assert cdr1as_sam.exists(), "cdr1as_anchors.sam missing"
    assert cdr1as_genome.exists(), "CDR1as_locus.fa missing"

    runner = CliRunner()
    result = runner.invoke(
        cli_main,
        [
            "call",
            str(cdr1as_sam),
            "--genome",
            str(cdr1as_genome),
            "--name",
            "cdr1as",
            "--prefix",
            "cdr1as_",
            "--anchor",
            "20",  # matches the typical anchor size used upstream
        ],
    )

    # Engine should run successfully
    assert result.exit_code == 0, result.output

    # Collect non-empty, non-comment lines from CLI output
    out_lines = [
        line.strip()
        for line in result.output.splitlines()
        if line.strip() and not line.startswith("#")
    ]

    # We only assert that at least one junction is marked CIRCULAR for this
    # CDR1as anchor test. Later, we can tighten this to check coordinates.
    circ_lines = []
    for ln in out_lines:
        fields = ln.split("\t")
        if not fields:
            continue
        categories = fields[-1].split(",")
        if "CIRCULAR" in categories:
            circ_lines.append(ln)

    assert circ_lines, "No CIRCULAR junctions emitted for CDR1as test"
