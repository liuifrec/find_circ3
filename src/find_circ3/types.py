# src/find_circ3/types.py
from __future__ import annotations

from typing import Protocol, runtime_checkable
import pysam


@runtime_checkable
class AlignedSegmentLike(Protocol):
    """
    Minimal interface we need from pysam.AlignedSegment.

    This is used purely for typing and duck-typing; anything that implements
    these attributes will be accepted by the breakpoint/accumulator code.
    """

    # Read / template id
    @property
    def query_name(self) -> str:
        ...

    # Genomic location
    @property
    def reference_name(self) -> str:
        ...

    @property
    def reference_start(self) -> int:
        ...

    @property
    def query_alignment_length(self) -> int:
        ...

    # Flags / quality
    @property
    def is_unmapped(self) -> bool:
        ...

    @property
    def is_reverse(self) -> bool:
        ...

    @property
    def mapping_quality(self) -> int:
        ...

    @property
    def cigarstring(self) -> str | None:
        ...


class GenomeAccessor:
    """
    Tiny adapter around pysam.FastaFile to match the expectations in
    breakpoints.py and engine.py.
    """

    def __init__(self, fasta_path: str) -> None:
        self._fa = pysam.FastaFile(fasta_path)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """
        Return reference sequence slice [start, end) for the given chromosome.
        """
        return self._fa.fetch(chrom, start, end)

    def close(self) -> None:
        self._fa.close()

    def __enter__(self) -> "GenomeAccessor":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()
