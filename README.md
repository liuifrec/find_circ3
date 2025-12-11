# find_circ3

A modern, Python 3 re‑implementation of Marvin Jens’ **find_circ.py v1.2** for circular RNA (circRNA) detection.

`find_circ3` keeps the original logic and output format for legacy anchor files, while providing a cleaner CLI and tests. It is designed to:

- Reproduce the original `find_circ.py` behavior on legacy unmapped‑anchor workflows.
- Work as a library/CLI component for downstream tools (e.g. `circyto`).
- Be fully testable and maintainable (Python 3, `pytest`, `pysam`).

> **Status (Dec 2025)**  
> - Legacy engine ported almost verbatim from `find_circ.py`.  
> - Regression tests against the original tiny test and CDR1as locus.  
> - Regression on a chr21 subset of ERR2139486: **53/54** original junctions recovered (the remaining one is a very long circular candidate spanning most of chr21).


---

## 1. Installation

Clone the repository and install in editable mode (recommended during development):

```bash
git clone https://github.com/liuifrec/find_circ3.git
cd find_circ3

# create & activate a virtualenv if you like
python -m venv .venv
source .venv/bin/activate

pip install -e .
```

This installs the CLI entry point:

```bash
find-circ3 --help
```

You should see the two subcommands:

- `find-circ3 anchors` – generate A/B anchors from **unmapped** reads.
- `find-circ3 call` – call circular and linear junctions from an **anchors SAM/BAM** aligned to the genome.


---

## 2. Concepts and input files

The pipeline follows the original find_circ workflow:

1. Map raw paired‑end reads to the genome.
2. Collect **unmapped reads** from the primary alignment.
3. Convert unmapped reads to **anchor reads** (left/right end “A/B” anchors).
4. Map anchors back to the genome.
5. From anchor‑anchor pairs, infer:
   - **Circular junctions** (back‑splice, reversed orientation).
   - **Linear junctions** (canonical splicing).

Key files:

- **Raw reads**: `sample_R1.fastq.gz`, `sample_R2.fastq.gz`
- **Genome FASTA**: `genome.fa` (and a bowtie2 index)
- **Primary alignment**: `sample_vs_genome.bam`
- **Unmapped reads BAM**: `sample_unmapped.bam`
- **Anchors FASTQ**: `sample_anchors.fastq.gz`
- **Anchors vs genome SAM/BAM**: `sample_anchors_vs_genome.sam`
- **Final junctions BED**: `sample_splice_sites.bed`


---

## 3. Full pipeline: FASTQ → circRNA junctions

Below is a minimal, end‑to‑end example using bowtie2 + samtools + find_circ3.  
Adjust paths, thread counts, and options as appropriate for your data.

### 3.1 Build bowtie2 index (once per reference)

```bash
bowtie2-build genome.fa genome_bt2
```

### 3.2 Map raw reads to the genome (primary mapping)

> ⚠️ **Important:** do **not** use `--no-unal` here – we need unmapped reads.

```bash
bowtie2 -p 8 --score-min C,-15,0 --mm --reorder   -x genome_bt2   -1 sample_R1.fastq.gz   -2 sample_R2.fastq.gz   -S sample_vs_genome.sam

# Convert to BAM and clean up
samtools view -bS sample_vs_genome.sam | samtools sort -o sample_vs_genome.bam
samtools index sample_vs_genome.bam
rm sample_vs_genome.sam
```

### 3.3 Extract unmapped reads

We keep only reads that **failed to align** in the primary mapping:

```bash
samtools view -b -f 4 sample_vs_genome.bam > sample_unmapped.bam
```

- `-f 4` → select unmapped reads  
- You can add more flags (e.g. to require both mates unmapped) if desired, but this simple form works well for the classic pipeline.

### 3.4 Generate anchor reads

Now convert unmapped reads into left/right anchors. The `anchors` command is a thin wrapper around the classic `unmapped2anchors` logic.

```bash
find-circ3 anchors sample_unmapped.bam   --anchor 20   --out sample_anchors.fastq.gz
```

Notes:

- `--anchor` must match the anchor size you will later pass to `find-circ3 call` (default 20 nt).
- Output is a FASTQ file where each original read yields a pair of anchors with `*_A` / `*_B` style names.
- Use `find-circ3 anchors --help` for additional options (e.g. compression, filtering).

### 3.5 Map anchors back to the genome

We now align the **anchors** to the genome. Here `--no-unal` is safe: unmapped anchors are simply ignored later.

```bash
bowtie2 -p 8 --score-min C,-15,0 --mm --reorder --no-unal   -x genome_bt2   -U sample_anchors.fastq.gz   -S sample_anchors_vs_genome.sam
```

You can optionally compress to BAM:

```bash
samtools view -bS sample_anchors_vs_genome.sam > sample_anchors_vs_genome.bam
```

Either SAM or BAM can be passed to `find-circ3 call`.

### 3.6 Call circular and linear junctions

This step reproduces the original `find_circ.py` behavior (for legacy anchor names) using a faithful Python 3 port under the hood.

```bash
find-circ3 call sample_anchors_vs_genome.sam   --genome genome.fa   --name sample   --prefix sample_   --anchor 20   --min-as-xs 2   --max-intron 200000   --min-support 1   > sample_splice_sites.bed
```

Key options:

- `--genome` / `-g` – reference FASTA (same as used for bowtie2 alignment).
- `--name` / `-n` – sample name used in the “tissues” column.
- `--prefix` / `-p` – prefix for junction IDs (`sample_circ_000001`, etc.).
- `--anchor` – **must match** anchor size used in `find-circ3 anchors`.
- `--min-as-xs` – minimal **AS–XS** margin to treat an anchor as uniquely placed (legacy `min_uniq_qual`, default 2).
- `--max-intron` – maximum genomic span allowed for a junction (default 200 kb).
- `--min-support` – minimal number of supporting reads per junction (default 1).
- `--allow-non-canonical` – if set, also considers non‑GT/AG motifs (not recommended for default runs).

Output (`sample_splice_sites.bed`) is a BED‑like table:

```text
chrom  start  end   name          n_reads strand n_uniq uniq_bridges best_qual_left best_qual_right tissues ...
chr21  33523435 33527305 sample_circ_000001 31 + 13 28 40 40 sample 31 ... CIRCULAR,PERFECT_EXT,ANCHOR_UNIQUE, ...
```

The **last column** lists comma‑separated categories, such as:

- `CIRCULAR` or `LINEAR`
- `CANONICAL` (GT/AG)
- `ANCHOR_UNIQUE`
- `PERFECT_EXT` / `GOOD_EXT` / `OK_EXT`
- `UNAMBIGUOUS_BP`

You can filter high‑confidence circRNAs by requiring, for example:

- `CIRCULAR`
- `ANCHOR_UNIQUE`
- `CANONICAL`
- `PERFECT_EXT` or `GOOD_EXT`


---

## 4. Tiny test and CDR1as regression checks

The repository ships small regression datasets that are also used by the test suite.

### 4.1 Tiny test

```bash
find-circ3 call tests/test_against_legacy/data/tiny.sam   --genome tests/test_against_legacy/data/tiny_genome.fa   --name test   --prefix test_   --anchor 5   > tiny_out.bed

diff <(grep -v '^#' tiny_out.bed)      <(grep -v '^#' tests/test_against_legacy/data/tiny_expected.bed)
```

This should produce identical non‑comment lines.

### 4.2 CDR1as locus

```bash
find-circ3 call tests/test_against_legacy/data/cdr1as_anchors.sam   --genome tests/test_against_legacy/data/cdr1as_genome.fa   --name cdr1as_test   --prefix cdr1as_   --anchor 20   > cdr1as_splice_sites.bed
```

You should see a clear circular signal at the CDR1as locus.

### 4.3 ERR2139486 chr21 subset (legacy regression)

In development, a chr21 subset of ERR2139486 was used to tune the legacy engine port. Running:

```bash
TEST=tests/test_data/ERR2139486_chr21

find-circ3 call   "$TEST/ERR2139486.anchors.sam"   --genome tests/data/ref/chr21.fa   --name ERR2139486_chr21   --prefix ERR2139486_   --anchor 20   > "$TEST/ERR2139486_splice_sites.bed"
```

produces **53/54** of the original junctions from the historical `find_circ.py` run on this subset. The missing junction is a very long circular candidate that spans most of chr21; everything else matches closely.


---

## 5. Engine behavior and compatibility notes

### 5.1 Two internal engines

`find_circ3 call` automatically chooses one of two internal engines:

1. **Legacy engine** (default for “real” anchors):
   - Triggered when QNAMEs look like `*_A__SEQ` / `*_B` (output of `find-circ3 anchors` or original `unmapped2anchors.py`).
   - Uses a faithful Python 3 port of `find_circ.py` (including breakpoint scanning, AS/XS‑based uniqueness, and category labels).
   - Used for real data and regression tests (CDR1as, ERR2139486 chr21).

2. **Modern HitAccumulator engine** (for simple SAM fixtures):
   - Triggered when input SAM/BAM has plain QNAMEs and no legacy anchor labels.
   - Used mainly for small synthetic tests like `tiny.sam`, where anchors are already “baked in”.

You generally do not need to care which engine is used – the choice is automatic based on QNAME patterns.

### 5.2 Things to avoid

- **Do not use `--no-unal` on the primary mapping.**  
  You *need* the unmapped reads to feed into `find-circ3 anchors`.

- **Do not feed primary alignment SAM/BAM directly into `find-circ3 call`.**  
  `find-circ3 call` expects **anchor alignments**, not arbitrary primary mappings.

- **Do not change the anchor size between steps.**  
  `--anchor` in `find-circ3 anchors` and `find-circ3 call` must match.


---

## 6. Development and tests

To run the tests (including legacy regressions):

```bash
pytest -q
```

The test suite covers:

- CLI regression on a tiny synthetic example.
- CLI regression on the CDR1as locus.
- Legacy breakpoint logic.
- Anchors CLI behavior.
- Legacy regression on ERR2139486 chr21 (large test, may be xfailed or omitted in CI).


---

## 7. License and attribution

- The original `find_circ.py` (c) Marvin Jens, 2012–2015.  
- `find_circ3` re‑implements that logic in Python 3, with additional tests and infrastructure for modern workflows.

Please cite the original **find_circ** paper and Marvin Jens’ work when using this tool in publications.
