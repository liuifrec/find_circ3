# src/find_circ3/cli.py
from __future__ import annotations

from pathlib import Path
from typing import Optional

import sys
import click

from .engine import run_find_circ
from .anchors import main as anchors_cmd  # click.command in anchors.py


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def main() -> None:
    """find_circ3 command-line interface.

    Subcommands:
      call    : run the core find_circ3 engine on an anchors SAM/BAM
      anchors : generate A/B anchors from unmapped reads
    """
    # Click uses this as the root group; no body needed.
    pass


@main.command(name="call")
@click.argument(
    "anchors_sam",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
)
@click.option(
    "--genome",
    "-g",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    help="Reference genome FASTA.",
)
@click.option(
    "--name",
    "-n",
    "sample_name",
    default="unknown",
    show_default=True,
    help="Sample name used in the 4th BED column.",
)
@click.option(
    "--prefix",
    "-p",
    default="",
    show_default=True,
    help="Prefix for junction names.",
)
@click.option(
    "--anchor",
    "-a",
    "anchor_size",
    default=20,
    show_default=True,
    help="Anchor size (must match unmapped2anchors3).",
)
@click.option(
    "--min-uniq-qual",
    default=2,
    show_default=True,
    help="Minimum AS–XS margin to treat an anchor as unique.",
)
@click.option(
    "--margin",
    type=int,
    default=None,
    help="Breakpoint search flank margin (default = anchor/4).",
)
@click.option(
    "--max-mismatches",
    default=2,
    show_default=True,
    help="Maximum allowed mismatches in breakpoint search.",
)
@click.option(
    "--strandpref/--no-strandpref",
    default=False,
    show_default=True,
    help="Prefer strand-matched breakpoints when tie-breaking.",
)
@click.option(
    "--stats",
    "stats_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional path to write run statistics.",
)
@click.option(
    "--reads",
    "reads_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional path to write supporting read sequences.",
)
def call_cmd(
    anchors_sam: Path,
    genome: Path,
    sample_name: str,
    prefix: str,
    anchor_size: int,
    min_uniq_qual: int,
    margin: Optional[int],
    max_mismatches: int,
    strandpref: bool,
    stats_path: Optional[Path],
    reads_path: Optional[Path],
) -> None:
    """Call circular / linear junctions from an anchors SAM/BAM file."""
    for line in run_find_circ(
        anchors_fastq=anchors_sam,
        genome=genome,
        sample_name=sample_name,
        prefix=prefix,
        min_uniq_qual=min_uniq_qual,
        anchor_size=anchor_size,
        stats_path=stats_path,
        reads_path=reads_path,
        margin=margin,
        max_mismatches=max_mismatches,
        strandpref=strandpref,
    ):
        click.echo(line)


# Add the anchors subcommand as an alias to the anchors.py CLI
# anchors_cmd is itself a @click.command() named "main" in anchors.py.
main.add_command(anchors_cmd, name="anchors")


if __name__ == "__main__":
    # This still works when running `python -m find_circ3.cli`
    sys.exit(main())
