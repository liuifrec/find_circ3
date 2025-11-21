# find_circ3 benchmarks

This document describes how to benchmark **find_circ3** against the original Python 2 **find_circ** implementation and how to reproduce the small regression tests used during development.

The idea is to keep a small, well‑documented collection of datasets and comparison scripts that can be cited in a methods section as:

> “find_circ3 faithfully reproduces the junction calls of the original find_circ v1.2 on the CDR1as locus, HEK293, and synthetic circRNA benchmarks.”

---

## 1. Repository layout (suggested)

In your `find_circ3` repo, you can use a layout like:

```text
benchmarks/
  legacy/
    README.md          # this document (or a variant of it)
    cdr1as/
      find_circ/       # outputs from legacy Python 2 find_circ
      find_circ3/      # outputs from find_circ3
    hek293/
      find_circ/
      find_circ3/
    synthetic/
      find_circ/
      find_circ3/
  reports/
    summary.tsv        # high‑level concordance metrics
```

You do not need all datasets initially; starting with **cdr1as** is already valuable.

---

## 2. Prerequisites

You need:

- A working Python 2 environment with the original `find_circ.py` and its helper scripts (`unmapped2anchors.py`, `cmp_bed.py`, etc.).
- A Python 3 environment with `find_circ3` installed (for example):

  ```bash
  cd find_circ3
  python -m venv .venv
  source .venv/bin/activate
  pip install -e .
  ```

- Aligners (as in the legacy test data), typically:
  - `bowtie2`
  - `samtools`

---

## 3. CDR1as parity benchmark

### 3.1 Run legacy find_circ (Python 2)

From the original `find_circ/test_data` directory:

```bash
# 1) build bowtie2 index (once)
bowtie2-build CDR1as_locus.fa bt2_cdr1as_locus

# 2) align reads
bowtie2 -p8 --very-sensitive --score-min=C,-15,0 --reorder --mm     -f -U cdr1as_reads.fa -x bt2_cdr1as_locus     2> bt2_firstpass.log   | samtools view -hbuS -   | samtools sort -o cdr1as_test.bam

# 3) fetch unmapped reads
samtools view -hf 4 cdr1as_test.bam | samtools view -Sb - > unmapped_cdr1as_test.bam

# 4) generate anchors (Python 2)
python2 ../unmapped2anchors.py unmapped_cdr1as_test.bam > anchors_cdr1as_test.fastq

# 5) align anchors and run legacy find_circ
bowtie2 -q -U anchors_cdr1as_test.fastq -x bt2_cdr1as_locus --reorder --mm --very-sensitive --score-min=C,-15,0     2> bt2_secondpass.log   | python2 ../find_circ.py -G CDR1as_locus.fa -n cdr1as_legacy -p cdr1as_         --stats legacy_sites.log         --reads legacy_spliced_reads.fa     > legacy_splice_sites.bed
```

Copy `legacy_splice_sites.bed` into:

```text
benchmarks/legacy/cdr1as/find_circ/splice_sites.bed
```

### 3.2 Run find_circ3 on the same anchors

Still in `find_circ3` repo (Python 3 venv active):

```bash
find-circ3 call   tests/test_against_legacy/data/cdr1as_anchors.sam   --genome tests/test_against_legacy/data/CDR1as_locus.fa   --name cdr1as_find_circ3   --prefix cdr1as_   --anchor 20   > benchmarks/legacy/cdr1as/find_circ3/splice_sites.bed
```

> If you prefer to regenerate `cdr1as_anchors.sam` yourself, you can repeat the bowtie2 alignment on `anchors_cdr1as_test.fastq` and save the SAM output instead of piping directly into `find_circ`.

### 3.3 Compare outputs

Use the original `cmp_bed.py` where possible, or simple `diff`/`awk` checks:

```bash
cd benchmarks/legacy/cdr1as

# direct text diff (strict)
diff -u find_circ/splice_sites.bed find_circ3/splice_sites.bed || true

# or compare only key columns (chrom, start, end, name, category)
awk '{print $1,$2,$3,$4,$18}' find_circ/splice_sites.bed > legacy_keycols.txt
awk '{print $1,$2,$3,$4,$18}' find_circ3/splice_sites.bed > fc3_keycols.txt
diff -u legacy_keycols.txt fc3_keycols.txt || true
```

You can collect concordance into a small `reports/summary.tsv` with columns like:

- dataset
- n_legacy_junctions
- n_find_circ3_junctions
- n_shared_exact
- n_shared_coordinate_only
- Jaccard_index

---

## 4. HEK293 and synthetic benchmarks

Once CDR1as parity looks good, you can repeat similar steps for:

### 4.1 HEK293 example

Using the legacy `hek_test` and `hek_test2` targets in `find_circ/test_data/Makefile`, capture:

```text
benchmarks/legacy/hek293/find_circ/splice_sites.bed
```

Then run:

```bash
find-circ3 call   <hek_anchors.sam or BAM converted to SAM>   --genome /path/to/hg19.fa   --name hek_find_circ3   --prefix hek_   --anchor 20   > benchmarks/legacy/hek293/find_circ3/splice_sites.bed
```

Finally, compare as in the CDR1as example.

### 4.2 Synthetic circRNA data

The legacy test data includes a synthetic read generator (e.g. from `wbcel235.circs.bed.gz`). If you choose to use it, you can add:

```text
benchmarks/legacy/synthetic/find_circ/splice_sites.bed
benchmarks/legacy/synthetic/find_circ3/splice_sites.bed
```

and compute the same concordance metrics.

---

## 5. Suggested text for a methods section

Once you are satisfied with parity, you can summarise the benchmark as:

> “find_circ3 is a Python 3 reimplementation of the original find_circ circRNA detector (Memczak et al., 2013). We confirmed that find_circ3 reproduces the junction calls of the original Python 2 implementation on the CDR1as locus, HEK293 test data, and synthetic circRNA benchmarks provided with the original package, matching coordinates and categories for >99% of junctions.”

You can tailor the exact wording once you have concrete numbers in `reports/summary.tsv`.
