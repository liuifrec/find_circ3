# find_circ3
*Modern Python 3 reimplementation of find_circ (2013)*

find_circ3 is a modern, actively maintained Python 3 reimplementation of the classic circRNA detector **find_circ** (Memczak et al., Nature 2013).

It preserves the algorithmic behavior of **find_circ v1.2** while providing:

- Python 3–native code  
- Clean modular API (compatible with **circyto**)  
- Unit-tested breakpoint search, edit-distance scoring, and anchor-overlap logic  
- Fully rewritten **unmapped2anchors.py** in Python 3  
- Deterministic, testable outputs  
- Simple CLI using *click*  
- Modern packaging via `pyproject.toml`  

---

## ✨ Features

### ✔ Faithful Reimplementation of find_circ
- Anchor pairing & back-splice geometry  
- True breakpoint detection (read-based splice reconstruction)  
- Edit distance + anchor overlap scoring  
- Canonical / noncanonical splice signal recognition  
  (GT/GC/AT–AG)  
- Strand-aware tie-breaking (optional)  
- AS/XS-based unique-bridge scoring  
- Full legacy category reconstruction:

```
CANONICAL, STRANDMATCH,
ANCHOR_UNIQUE, NO_UNIQ_BRIDGES,
UNAMBIGUOUS_BP,
PERFECT_EXT, GOOD_EXT, OK_EXT,
CIRCULAR or LINEAR
```

---

## ✔ Modular Architecture

```text
src/find_circ3/
  cli.py             # 'find-circ3' CLI (group: call, anchors)
  anchors.py         # Python 3 unmapped2anchors3 (find-circ3-anchors / find-circ3 anchors)
  engine.py          # main junction detection logic
  breakpoints.py     # breakpoint detection + scoring
  hit_accumulator.py # scoring + category assembly
  io.py              # alignment helpers
  version.py
```

---

## ✔ High-Fidelity Tests

Includes reproducible tests for:

- tiny SAM → tiny expected BED  
- CDR1as anchor detection  
- Breakpoint detection (edit distance, overlap)  
- Anchors CLI tests (BAM + FASTQ modes)  
- Legacy dataset consistency  
- Optional integration test: **BAM → anchors → bowtie2 → call**

All tests pass in the latest repo state.

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

## 🧬 End-to-End Workflow  
**Paired-end FASTQ → circRNA BED**  
*(Fully working, no placeholders)*

This reproduces the original find_circ pipeline entirely in Python 3.

---

### Step 0 — Inputs

**Paired-end FASTQ**
```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

**Reference genome**
```text
genome.fa
```

**Bowtie2 index**
```text
genome_index
```

---

### Step 1 — Map paired-end reads & extract unmapped

Unmapped reads become anchor candidates.  
**Important:** do *not* use `--no-unal` here, we need the unmapped reads in the BAM.

```bash
bowtie2 -x genome_index   -1 sample_R1.fastq.gz   -2 sample_R2.fastq.gz   --very-sensitive   --score-min C,-15,0   2> sample_firstpass.log   | samtools view -bS - > sample_aln.bam

samtools view -b -f 4 sample_aln.bam > sample_unmapped.bam
```

---

### Step 2 — Generate anchors (Python 3 unmapped2anchors3)

You can now use **either** the dedicated console script **or** the grouped CLI alias:

```bash
# Recommended: grouped CLI
find-circ3 anchors sample_unmapped.bam   --anchor 20   --min-qual 5   > sample_anchors.fastq

# Exact equivalent legacy-style entry point
find-circ3-anchors sample_unmapped.bam   --anchor 20   --min-qual 5   > sample_anchors.fastq
```

**Output example:**

```text
@READ123_A__ACCGTACTGAGT...
ACCGTACTGAGT...
+
IIIIIIIIIIIIIII

@READ123_B
...TGATCGTTACGA
+
IIIIIIIIIIIIIII
```

- **A-anchor** embeds full read sequence  
- **B-anchor** references the same read without embedding  
- Exactly matches original find_circ format  

---

### Step 3 — Map anchors

```bash
bowtie2 -q   -U sample_anchors.fastq   -x genome_index   --reorder --mm --very-sensitive   --score-min C,-15,0   2> sample_secondpass.log   > sample_anchors.sam
```

Output must include:  
- proper SAM alignments  
- **AS** and **XS** tags (Bowtie2 default)  

---

### Step 4 — Call circRNA junctions

```bash
find-circ3 call sample_anchors.sam   --genome genome.fa   --name sample1   --prefix sample1_   --anchor 20   --min-uniq-qual 2   --max-mismatches 2   --margin 5   --strandpref   --stats sample_sites.log   --reads sample_spliced_reads.fa   > sample_splice_sites.bed
```

Key options:

| Option | Meaning |
|--------|---------|
| `--anchor`          | Must match Step 2                      |
| `--min-uniq-qual`  | AS–XS threshold                        |
| `--max-mismatches` | Allowed mismatches                     |
| `--margin`         | Flank margin (default = anchor/4)      |
| `--strandpref`     | Tie-breaking by strand                 |
| `--stats`          | Write stats log                        |
| `--reads`          | Write supporting read sequences        |

---

## 🧾 Output Format (18-column BED-like)

Example:

```text
chr1  5  22  sample1_chr1:5|22  2  +  2  1  60  60  sample  2  0  0  1  GTAG  MATCH  ANCHOR_UNIQUE,CANONICAL,CIRCULAR,PERFECT_EXT,STRANDMATCH,UNAMBIGUOUS_BP
```

Columns:

1. chrom  
2. start  
3. end  
4. name  
5. n_reads  
6. strand  
7. n_uniq  
8. uniq_bridges  
9. best_qual_left  
10. best_qual_right  
11. tissues  
12. tiss_counts  
13. edits  
14. anchor_overlap  
15. breakpoints  
16. signal  
17. strandmatch  
18. category list  

Output is deterministic and unit-tested.

---

## 🧪 Bundled Examples (Tests)

- `tiny.sam` → expected BED  
- CDR1as → produces ≥1 CIRCULAR call  
- breakpoints tests  
- anchors tests  
- optional **integration test**:
  - builds a tiny genome index,
  - creates a small unmapped BAM,
  - runs `find-circ3 anchors` → bowtie2 → `find-circ3 call`,
  - asserts that at least one junction line is produced.

---

## 🤝 Attribution

Original algorithm & Python 2 version:  
**Marvin Jens & Rajewsky Lab (2013)**

Python 3 modernization and reimplementation:  
**Yu-Chen (James) Liu**

---

## 📌 License

**GPL-3.0-or-later**  
(Fully compatible with the original licensing)

---

## 🧭 Roadmap

- ✔ Full breakpoint engine  
- ✔ AS/XS scoring + category reconstruction  
- ✔ Python 3 unmapped2anchors3  
- ✔ Consistent CLI (`--margin`, `--max-mismatches`, `--strandpref`)  
- ✔ `find-circ3 anchors` CLI alias  
- 🔄 circyto integration  
- 🔄 Multi-threading & performance optimization  
- 🔄 PyPI distribution  

---

## 💪 Status

**find_circ3 is now complete and stable:**

- Full anchor → alignment → junction pipeline  
- All find_circ v1.2 behaviors reproduced  
- No placeholders  
- Full test suite coverage  
