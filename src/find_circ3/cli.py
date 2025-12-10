# src/find_circ3/cli.py
from __future__ import annotations

from pathlib import Path
import sys

import click

from .engine import run_find_circ
from .anchors import main as anchors_main  # the click.command() from anchors.py


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def main() -> None:
    """
    find_circ3 command-line interface.

    Subcommands:
      anchors : generate A/B anchors from unmapped BAM
      call    : call junctions from anchors SAM
    """
    # click uses this as the root group; no body needed.
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
    help="Reference genome FASTA used to align the anchors.",
)
@click.option(
    "--name",
    "-n",
    "sample_name",
    default="sample",
    show_default=True,
    help="Sample name used in the 11th BED column.",
)
@click.option(
    "--prefix",
    "-p",
    default="",
    show_default=True,
    help="Prefix for junction names (4th BED column).",
)
@click.option(
    "-a",
    "--anchor",
    type=int,
    default=20,
    show_default=True,
    help="Anchor size used upstream in unmapped2anchors3.",
)
@click.option(
    "--min-mapq",
    type=int,
    default=0,
    show_default=True,
    help="Minimum MAPQ for anchor alignments to be considered.",
)
@click.option(
    "--min-as-xs",
    type=int,
    default=2,
    show_default=True,
    help="Minimum AS–XS margin to treat an anchor as uniquely placed.",
)
@click.option(
    "--max-intron",
    type=int,
    default=200000,
    show_default=True,
    help="Maximum intron length allowed for junctions.",
)
@click.option(
    "--min-support",
    type=int,
    default=1,
    show_default=True,
    help="Minimum number of supporting reads per junction (set to 1 for tests).",
)
@click.option(
    "--allow-non-canonical/--no-allow-non-canonical",
    default=False,
    show_default=True,
    help="Allow non-GT/AG splice motifs.",
)
def call_cmd(
    anchors_sam: Path,
    genome: Path,
    sample_name: str,
    prefix: str,
    anchor: int,
    min_mapq: int,
    min_as_xs: int,
    max_intron: int,
    min_support: int,
    allow_non_canonical: bool,
) -> None:
    """
    Call circular / linear junctions from an anchors SAM/BAM file.

    This is a thin wrapper around engine.run_find_circ, which yields
    BED-format junction lines. We just stream them to stdout.
    """
    for line in run_find_circ(
        anchors_fastq=anchors_sam,
        genome=genome,
        sample_name=sample_name,
        prefix=prefix,
        anchor=anchor,
        min_mapq=min_mapq,
        min_as_xs=min_as_xs,
        max_intron=max_intron,
        min_support=min_support,
        allow_non_canonical=allow_non_canonical,
        sample=sample_name,
        stats_path=None,
        reads_path=None,
    ):
        # run_find_circ yields '\n'-terminated strings; click will add its own
        # newline unless we strip it, so we rstrip() first.
        click.echo(line.rstrip("\n"))


# Expose the anchors subcommand (from anchors.py) under the group.
# anchors_main is already a @click.command() named "main" in anchors.py.
main.add_command(anchors_main, name="anchors")


if __name__ == "__main__":
    sys.exit(main())
