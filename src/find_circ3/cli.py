import sys
from pathlib import Path

import click

from .version import __version__
from .engine import run_find_circ


@click.group()
@click.version_option(__version__, prog_name="find-circ3")
def main():
    """find-circ3: Python 3 reimplementation of find_circ."""
    pass


@main.command("call")
@click.argument(
    "anchors_sam",
    type=click.Path(exists=True, readable=True),
)
@click.option(
    "--genome",
    "-G",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Reference genome FASTA (single multi-FASTA file).",
)
@click.option(
    "--name",
    "-n",
    default="unknown",
    show_default=True,
    help="Sample name / tissue name.",
)
@click.option(
    "--prefix",
    "-p",
    default="",
    show_default=True,
    help="Prefix to prepend to each junction name.",
)
@click.option(
    "--min-uniq-qual",
    "-q",
    default=2,
    show_default=True,
    help="Minimal uniqueness score for anchor alignments (AS–XS margin).",
)
@click.option(
    "--anchor",
    "-a",
    default=20,
    show_default=True,
    help="Anchor size (must match unmapped2anchors3).",
)
@click.option(
    "--stats",
    type=click.Path(writable=True),
    help="Write numeric run statistics to this file.",
)
@click.option(
    "--reads",
    type=click.Path(writable=True),
    help="Write supporting reads to this file instead of stderr.",
)
@click.option(
    "--margin",
    type=int,
    help="Breakpoint flank margin (bp); defaults to anchor_size // 4.",
)
@click.option(
    "--max-mismatches",
    "max_mismatches",
    type=int,
    default=2,
    show_default=True,
    help="Maximum mismatches allowed in breakpoint search.",
)
@click.option(
    "--strandpref/--no-strandpref",
    default=False,
    show_default=True,
    help="Prefer strand-matched breakpoints when breaking ties.",
)
def call_cmd(
    anchors_sam,
    genome,
    name,
    prefix,
    min_uniq_qual,
    anchor,
    stats,
    reads,
    margin,
    max_mismatches,
    strandpref,
):
    """
    Detect linear and circular junctions from anchor alignments.

    ANCHORS_SAM:
        SAM/BAM file of anchor alignments (e.g. bowtie2 output on anchors).
    """
    anchors_path = Path(anchors_sam)
    genome_path = Path(genome)
    stats_path = Path(stats) if stats else None
    reads_path = Path(reads) if reads else None

    try:
        junctions = run_find_circ(
            anchors_fastq=anchors_path,
            genome=genome_path,
            sample_name=name,
            prefix=prefix,
            min_uniq_qual=min_uniq_qual,
            anchor_size=anchor,
            stats_path=stats_path,
            reads_path=reads_path,
            margin=margin,
            max_mismatches=max_mismatches,
            strandpref=strandpref,
        )
        for bed_line in junctions:
            click.echo(bed_line)
    except Exception as e:
        # Generic safety net: report error and exit non-zero
        click.echo(f"# find-circ3 error: {e}", err=True)
        sys.exit(1)
