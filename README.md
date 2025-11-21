# find_circ3

A Python 3–compatible reimplementation and modernization of the classic
circRNA detector **find_circ** (Memczak et al., Nature 2013).

- Original algorithm and Python 2 implementation: Marvin Jens, Rajewsky Lab  
- Python 3 port, modernization, packaging, and circyto integration:
  **Yu-Chen (James) Liu**

## Status

**Phase 0** – packaging & CLI skeleton.  
Core detection engine is **not** implemented yet.

## Installation (development)

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scriptsctivate
pip install -e .
```

## Command-line interface

```bash
find-circ3 --help
find-circ3-anchors --help
```

These currently expose a basic CLI skeleton and stubbed functionality;
the core algorithm will be added in later phases.
