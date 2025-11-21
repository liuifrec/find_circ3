import sys
from pathlib import Path
import click
from .version import __version__

@click.group()
@click.version_option(__version__, prog_name="find-circ3")
def main():
    """find-circ3: Python 3 reimplementation of find_circ."""
    pass


@main.command("call")
@click.argument("anchors_fastq", type=click.Path(exists=True, readable=True))
@click.option("--genome", "-G", required=True, type=click.Path(exists=True, readable=True))
@click.option("--name", "-n", default="unknown")
@click.option("--prefix", "-p", default="")
@click.option("--min-uniq-qual", "-q", default=2)
@click.option("--anchor", "-a", default=20)
@click.option("--stats", type=click.Path(writable=True))
@click.option("--reads", type=click.Path(writable=True))
def call_cmd(anchors_fastq, genome, name, prefix, min_uniq_qual, anchor, stats, reads):
    """Detect junctions from anchor reads (stub)."""
    click.echo("# find-circ3 skeleton: engine not implemented yet", err=True)
    sys.exit(1)
