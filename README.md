# find_circ3

**find_circ3** is a modern, Python 3, actively maintained reimplementation of the classic circRNA detector **find_circ** (Memczak et al., *Nature* 2013).

It preserves the original algorithmic behavior of find_circ v1.2 while providing:

- Fully Python 3–native code  
- A clean modular API (for integration in pipelines like *circyto*)  
- Unit-tested breakpoint search, edit-distance scoring, and anchor-overlap logic  
- CLI built on `click`  
- Easier installation and packaging with `pyproject.toml`  
- Deterministic, testable outputs

---

## ✨ Features

### ✔ Faithful reimplementation of find_circ

- Anchor pairing and back-splice geometry  
- Real breakpoint detection (read-based splice reconstruction)  
- Edit distance, anchor overlap, GT/GC/AT–AG recognition  
- Canonical / noncanonical signal classification  
- Strand-aware tie-breaking (optional)  
- AS/XS–based unique-bridge scoring  
- Legacy category reconstruction:
  - `CANONICAL`, `STRANDMATCH`  
  - `ANCHOR_UNIQUE`, `NO_UNIQ_BRIDGES`  
  - `UNAMBIGUOUS_BP`  
  - `PERFECT_EXT`, `GOOD_EXT`, `OK_EXT`  
  - final `CIRCULAR` / `LINEAR` tag

### ✔ Modular architecture

```text
src/find_circ3/
  cli.py             # 'find-circ3' entry point
  anchors.py         # unmapped2anchors3 (Python 3 port)
  engine.py          # main junction detection logic
  breakpoints.py     # breakpoint search and scoring
  hit_accumulator.py # scoring + category assignment
  io.py              # BAM/SAM helpers (extendable)
  version.py
```

### ✔ High-fidelity testing

All logic is verified by a set of regression tests:

- **tiny.sam** minimal example  
- **CDR1as** landmark circRNA test  
- **legacy anchors** existence & non-emptiness  
- Dedicated **breakpoints** unit test validating edit distance + overlap  
- **anchors CLI** tests for BAM and FASTQ inputs

---

## 🚀 Installation

```bash
git clone https://github.com/liuifrec/find_circ3
cd find_circ3
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Run tests:

```bash
pytest -q
```

---

## 🧬 Typical workflow: from paired-end FASTQ to circRNA list

This section describes a full, reproducible workflow starting from **paired-end FASTQ** files and ending with a **circRNA junction BED** file, closely matching the original find_circ pipeline.

### Step 0 — Inputs

- Paired-end FASTQ files:  
  - `sample_R1.fastq.gz`  
  - `sample_R2.fastq.gz`
- Reference genome (FASTA):  
  - `genome.fa`
- Bowtie2 index built on the same reference:  
  - prefix: `genome_index`

> Adjust names and paths as needed.

---

### Step 1 — Map paired-end reads and collect unmapped reads

First, align the paired-end reads to the genome and keep **only unmapped reads**, which are the substrate for anchor generation.

```bash
# 1. Align paired-end reads
bowtie2 -x genome_index   -1 sample_R1.fastq.gz   -2 sample_R2.fastq.gz   --very-sensitive   --score-min C,-15,0   --no-unal   2> sample_firstpass.log   | samtools view -bS - > sample_aln.bam

# 2. Extract unmapped reads
samtools view -b -f 4 sample_aln.bam > sample_unmapped.bam
```

You can tune Bowtie2 options as you like; the above mirrors the original find_circ examples.

---

### Step 2 — Generate anchors (unmapped2anchors3)

Use the **find-circ3-anchors** command (Python 3 port of `unmapped2anchors.py`) to split each unmapped read into A/B anchors in FASTQ format.

```bash
find-circ3-anchors sample_unmapped.bam   --anchor 20   --min-qual 5   > sample_anchors.fastq
```

Key options:

- `--anchor` (`-a`): anchor size (typically 20–25).  
- `--min-qual` (`-q`): minimum average base quality for both anchors (legacy default is 5).  

The output `sample_anchors.fastq` contains reads like:

- `@READNAME_A__FULLSEQ` → left anchor  
- `@READNAME_B` → right anchor  

The full original read sequence is embedded in the A-read name (used for breakpoint reconstruction).

---

### Step 3 — Align anchors to the genome

Align the anchor FASTQ to the same reference genome and produce a SAM file. This mirrors the **second bowtie2 pass** in the legacy pipeline.

```bash
bowtie2 -q   -U sample_anchors.fastq   -x genome_index   --reorder --mm --very-sensitive   --score-min C,-15,0   2> sample_secondpass.log   > sample_anchors.sam
```

Notes:

- Output is **SAM on stdout**, written to `sample_anchors.sam`.  
- You can instead pipe SAM directly into `find-circ3` (see below) if you prefer a streamed workflow.

---

### Step 4 — Call circRNA junctions with find_circ3

Now run `find-circ3 call` on the **anchor alignment SAM** and the reference genome:

```bash
find-circ3 call sample_anchors.sam   --genome genome.fa   --name sample1   --prefix sample1_   --anchor 20   --min-uniq-qual 2   --max-mismatches 2   --margin 5   --strandpref   > sample_splice_sites.bed
```

Important options:

- `sample_anchors.sam` — SAM with anchor alignments (from Step 3).  
- `--genome` / `-G` — reference genome FASTA (`genome.fa`).  
- `--name` / `-n` — sample name used in stats and some fields.  
- `--prefix` / `-p` — prefix for junction names (e.g. `sample1_chr1:...`).  
- `--anchor` / `-a` — anchor size (must match Step 2).  
- `--min-uniq-qual` / `-q` — minimum AS–XS margin for anchors to count as “unique”.  
- `--margin` — breakpoint flank margin; defaults to `anchor_size // 4` if omitted.  
- `--max-mismatches` — maximum mismatches allowed in breakpoint search.  
- `--strandpref/--no-strandpref` — prefer strand-matched breakpoints when breaking ties.  

Optional outputs:

- `--stats sites.log` — write numeric summary statistics.  
- `--reads spliced_reads.fa` — write supporting reads used for splice site detection.

---

### One-line streamed variant (no intermediate SAM)

You can also stream Bowtie2 output directly into `find-circ3 call`:

```bash
bowtie2 -q   -U sample_anchors.fastq   -x genome_index   --reorder --mm --very-sensitive   --score-min C,-15,0   2> sample_secondpass.log   | find-circ3 call /dev/stdin         --genome genome.fa         --name sample1         --prefix sample1_         --anchor 20
```

Here `/dev/stdin` is used as the SAM path; this works on Unix-like systems.

---

## 📦 Output format

The output is a BED-like tab-separated file with 18 columns, matching the legacy `splice_sites.bed` format:

```text
1: chrom          (e.g. chr1)
2: start          (0-based)
3: end            (1-based)
4: name           (e.g. sample1_chr1:123|456)
5: n_reads        (# supporting reads)
6: strand         (+ or -)
7: n_uniq         (# anchors with unique placement)
8: uniq_bridges   (# reads where both anchors are unique)
9: best_qual_left
10: best_qual_right
11: tissues       (comma-separated; simplified to 'sample' by default)
12: tiss_counts   (comma-separated counts per tissue)
13: edits         (min edit distance at breakpoint)
14: anchor_overlap
15: breakpoints   (# tied best breakpoints)
16: signal        (e.g. GTAG)
17: strandmatch   (MATCH / MISMATCH / NA)
18: category      (comma-separated legacy-style labels)
```

Example line:

```text
chr1    5   22  sample1_chr1:5|22  2   +   2   1   60  60  sample  2  0  0  1  GTAG  MATCH  ANCHOR_UNIQUE,CANONICAL,CIRCULAR,PERFECT_EXT,STRANDMATCH,UNAMBIGUOUS_BP
```

---

## 🧪 Small example with bundled tests

The test suite ships with tiny synthetic data and a CDR1as example, demonstrating both the CLI and the breakpoint engine:

- `tests/test_against_legacy/test_regression.py`
  - runs `find-circ3 call` on a tiny SAM  
  - runs `find-circ3 call` on `cdr1as_anchors.sam` + `CDR1as_locus.fa`
- `tests/test_breakpoints.py`
  - validates edit-distance and overlap logic in isolation
- `tests/test_anchors_cli.py`
  - checks `find-circ3-anchors` behaviour for BAM and FASTQ inputs

You can use these as templates to design your own integration tests.

---

## 📁 Benchmarks

You can add a benchmark README under `benchmarks/` describing how to reproduce parity with the original Python 2 **find_circ** on:

- CDR1as locus  
- HEK293 example  
- Synthetic circRNA datasets  

These benchmarks can be summarized in a simple table and cited in a methods section.

---

## 🤝 Attribution

Original algorithm and Python 2 implementation:  
**Marvin Jens & Rajewsky Lab (2013).**

Python 3 reimplementation, modernization, packaging, and testing:  
**Yu-Chen (James) Liu**

---

## 📌 License

**GPL-3.0-or-later**, in full compliance with the original find_circ licensing.

---

## 🧭 Roadmap

- [x] Full breakpoint engine  
- [x] AS/XS scoring + category reconstruction  
- [x] `unmapped2anchors3` rewrite (Python 3)  
- [ ] Integration into *circyto* as a pluggable detector  
- [ ] Multi-threaded performance enhancements  
- [ ] Add conda & PyPI distribution  

---

## 💪 Notes

This tool is designed for scientific reproducibility and integration into modern pipelines:

- Deterministic, unit-tested core  
- Explicit, readable configuration via CLI or API  
- Clean separation between anchor generation, mapping, and junction calling  

Ideal for use as a detector backend in `circyto` and other circRNA analysis frameworks.
