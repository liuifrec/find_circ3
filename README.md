# find_circ3

find_circ3 is a fully rewritten high-accuracy circRNA detector for short‑read RNA‑seq, rebuilt with a modern breakpoint engine, AS/XS scoring, and a clean unmapped→anchors→calls workflow.

This README documents the **complete paired-end FASTQ → circRNA** workflow and all internal components.

---

## 🔧 Overview

**find_circ3** takes paired-end RNA‑seq FASTQs and identifies back-splice junctions (BSJs) using:

1. **unmapped2anchors3**  
   Extracts anchor candidates from the unmapped reads of a BWA-MEM alignment.
2. **find-circ3-anchors**  
   Aligns anchors, scores them (AS/XS), identifies candidate breakpoints.
3. **find-circ3-call**  
   Filters candidate breakpoints, enforces strand rules, reports final circRNAs.

All three components are now fully updated and pass tests.

---

## 📦 Pipeline Summary

### Step 1 — Align reads with BWA-MEM (retain unmapped)
```bash
bwa mem -T 19 -t 8 reference.fa reads_1.fq.gz reads_2.fq.gz > aln.sam
```

### Step 2 — Extract anchor sequences
```bash
find-circ3 unmapped2anchors3 aln.sam anchors.fa
```

### Step 3 — Align anchors to genome  
(uses BWA‑MEM internally and produces *.anchors.sam*)
```bash
find-circ3 anchors anchors.fa > anchors.sam
```

### Step 4 — Call circRNAs from aligned anchors
```bash
find-circ3 call anchors.sam circ_calls.txt
```

---

## 📁 Output Files

| File | Description |
|------|-------------|
| `anchors.fa` | 20–25 bp anchor sequences extracted from unmapped reads |
| `anchors.sam` | BWA‑MEM alignments of anchors back to genome |
| `circ_calls.txt` | Final circRNA calls after scoring/filtering |

---

## 🧠 Detection Logic

### 1. Anchor extraction (unmapped2anchors3)
- Reads unmapped fragments from SAM.
- Extracts both ends of each read (configurable anchor length).
- Writes anchors as FASTA.

### 2. Anchor alignment (find-circ3-anchors)
Each anchor receives:
- **AS score** (match quality)
- **XS score** (suboptimal match)
- **XS/AS ratio** filter  
- Left/right anchor pairing logic
- Preliminary breakpoint inference

### 3. Final BSJ calling (find-circ3-call)
Enforces:
- Same‑chromosome constraint  
- Valid orientation for a back-splice  
- Genomic distance limits  
- AS/XS score thresholds  
- Deduplication of identical junctions  

---

## 🐍 Installation

Direct executable via Python:

```bash
pip install -r requirements.txt
python setup.py install
```

Or as CLI:

```bash
find-circ3 --help
```

---

## 🧪 Tests

All tests under `tests/` cover:
- unmapped2anchors3 rewrite
- anchor scoring & XS/AS logic
- breakpoint inference
- full FASTQ→SAM→anchors→call workflow

Run tests:

```bash
pytest -q
```

---

## 🗺 Roadmap for circyto Integration

- Add detector API wrapper  
- Create consistent circ_feature_table.tsv  
- Provide find_circ3 outputs for multimodal export  
- Benchmark vs CIRI‑full / find_circ / CIRCexplorer2

---

## ✨ Citation

If you use **find_circ3**, please cite:

> Liu Y.-C., find_circ3 project (2025). High‑accuracy short‑read circRNA detection via breakpoint scoring and anchor‑based reconstruction.

---

