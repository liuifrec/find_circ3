import click
import pysam

@click.command()
@click.argument("input_bam", type=click.Path(exists=True, readable=True))
def main(input_bam):
    """Stub for unmapped2anchors."""
    bam = pysam.AlignmentFile(input_bam, "rb")
    unmapped = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            unmapped += 1
    bam.close()
    click.echo(f"# find-circ3 anchors stub: {unmapped} unmapped reads", err=True)


if __name__ == "__main__":
    main()
