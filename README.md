# find_circ3

**find_circ3** is a modern, Python 3, actively maintained reimplementation of the classic circRNA detector **find_circ** (Memczak et al., *Nature* 2013).

It preserves the original algorithmic behavior of find_circ v1.2 while providing:

- Fully Python 3–native code  
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
  - `CANONICAL`, `STRANDMATCH`,  
  - `ANCHOR_UNIQUE`, `NO_UNIQ_BRIDGES`,  
  - `UNAMBIGUOUS_BP`,  
  - `PERFECT_EXT`, `GOOD_EXT`, `OK_EXT`,  
  - and final `CIRCULAR` / `LINEAR`

### ✔ Fully modular architecture
```
find_circ3/
  engine.py          # main junction detection logic
  breakpoints.py     # full breakpoint search
  hit_accumulator.py # scoring + category assignment
  anchors.py         # future unmapped2anchors3 interface
  cli.py             # 'find-circ3' command
  io.py              # BAM/SAM helpers (expandable)
```

### ✔ High-fidelity testing  
All logic is verified by a set of regression tests:

- **tiny.sam** minimal example  
- **CDR1as** landmark circRNA test  
- **legacy files exist + non-empty**  
- A dedicated **breakpoints** unit test validating real edit distance + overlap

---

## 🚀 Installation

```bash
git clone https://github.com/liuifrec/find_circ3
cd find_circ3
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Test installation:

```bash
pytest -q
```

---

## 🧬 Usage

### Basic circular junction calling:

```bash
find-circ3 call anchors.sam   --genome genome.fa   --name sample1   --prefix samp1_   --anchor 20
```

### Options

| Option | Meaning |
|--------|----------|
| `--genome` | Reference genome FASTA |
| `--anchor` | Anchor size (default 20) |
| `--min-uniq-qual` | Minimum AS–XS margin for unique anchors |
| `--margin` | Breakpoint search flank margin (default anchor/4) |
| `--max-mismatches` | Max allowed mismatches in breakpoint search |
| `--strandpref/--no-strandpref` | Apply strand-aware tie-breaking |
| `--stats` | Write stats file |
| `--reads` | Write supporting read sequences |

---

## 📦 Example output (BED-like)

The output matches the legacy 18-column `splice_sites.bed` format, including composite categories:

```
chr1    5   22  samp1_chr1:5|22   2   +   2   1   60  60  sample  2  0  0  1  GTAG  MATCH  CANONICAL,PERFECT_EXT,STRANDMATCH,ANCHOR_UNIQUE,UNAMBIGUOUS_BP,CIRCULAR
```

---

## 📁 Benchmarks

A benchmark README is provided in:

```
benchmarks/find_circ3_benchmarks_README.md
```

It describes how to reproduce parity with original find_circ on:

- CDR1as locus  
- HEK293 example  
- Synthetic circRNA datasets  

---

## 🤝 Attribution

Original algorithm and Python 2 implementation:  
**Marvin Jens & Rajewsky Lab (2013).**

Python 3 reimplementation, modernization, packaging, and testing:  
**Yu-Chen (James) Liu**

---

## 📌 License

**GPL-3.0-or-later**, in full compliance with original find_circ licensing.

---

## 🧭 Roadmap

- [x] Full breakpoint engine  
- [x] AS/XS scoring + category reconstruction  
- [ ] Integration into *circyto* as a pluggable detector  
- [ ] `unmapped2anchors3` rewrite  
- [ ] Multi-threaded performance enhancements  
- [ ] Add conda & PyPI distribution  

---

## 💪 Notes

This tool was built with scientific reproducibility in mind.  
Every component is unit-tested, deterministic, and ready for integration in modern circRNA analysis workflows.

