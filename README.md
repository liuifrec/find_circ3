
# find_circ3

A modern, Python 3 reimplementation of **find_circ.py** with:

- a small, testable codebase,
- a dual-engine design (strict legacy mode + modern engine),
- and regression tests that lock behavior to the original `find_circ.py` and CDR1as examples.

The goal is **practical, stable circRNA calling** that behaves like the original tool on legacy anchor files, while still being hackable and maintainable.

---

## 1. Installation

```bash
git clone https://github.com/liuifrec/find_circ3.git
cd find_circ3

# Create a fresh virtualenv (recommended)
python3 -m venv .venv
source .venv/bin/activate

# Install in editable mode with all dependencies
pip install -U pip
pip install -e .
```

Requirements (roughly):

- Python ≥ 3.10
- `pysam`, `numpy`, `click`, `pytest`
- A short-read aligner (e.g. **bowtie2**)
- `samtools` for typical BAM/SAM processing

Run tests to confirm everything is wired correctly:

```bash
pytest -q
```

You should see:

- 6 tests **passed**
- 1 test **xpassed** (integration test that is intentionally skipped/expected)

---

## 2. Command line usage

The CLI exposes two subcommands:

- `find-circ3 anchors` – generate A/B anchors from unmapped reads  
- `find-circ3 call` – call circular and linear junctions from anchor alignments

### 2.1. `find-circ3 call`

This is the main entry point when you already have **anchor alignments** (SAM/BAM).

```bash
find-circ3 call \
  anchors.sam \
  --genome genome.fa \
  --name SAMPLE_NAME \
  --prefix SAMPLE_ \
  --anchor 20 \
  --max-intron 200000 \
  --min-support 1 \
  > splice_sites.bed
```

Key options:

- `anchors.sam`  
  SAM (or BAM) file produced by aligning A/B anchors against the genome.

- `--genome` / `-g`  
  Reference genome FASTA used for anchor alignment.

- `--name` / `-n`  
  Sample name written into the BED “tissues” column.

- `--prefix` / `-p`  
  Prefix for junction names, e.g. `SAMPLE_circ_000001`.

- `--anchor` / `-a`  
  Anchor size used upstream in the anchor-generation step  
  (must match the anchor length in `unmapped2anchors3` / `find-circ3 anchors`).

- `--max-intron`  
  Maximum allowed genomic span (`end - start`) in bp.  
  Default: **200,000**.  
  This is a *safety* filter to remove implausibly huge junctions on real data.

- `--min-support`  
  Minimum number of supporting reads per junction (default = 1).

- `--allow-non-canonical`  
  If set, also allow non-GT/AG splice motifs.  
  By default, only canonical **GT/AG** junctions are reported.

Run `find-circ3 call --help` to see the full set of options.

---

### 2.2. `find-circ3 anchors` (A/B anchor generation)

The `anchors` subcommand is the modern replacement for the original `unmapped2anchors.py` / `unmapped2anchors3` scripts:

```bash
find-circ3 anchors \
  unmapped.bam \
  --anchor 20 \
  > anchors.fastq
```

Typical pattern:

1. Start from **unmapped reads** from your original genome alignment.
2. Convert unmapped reads into A/B anchors (`anchors.fastq`).
3. Align these anchors to the genome with a sensitive local aligner (e.g. bowtie2).
4. Feed the resulting `anchors.sam` into `find-circ3 call`.

For the exact options and recommended parameters, always check:

```bash
find-circ3 anchors --help
```

> **Note:**  
> Earlier README versions accidentally suggested using `--no-unal` in anchor-alignment commands, which can silently strip out useful anchors. The current docs purposely avoid `--no-unal` here; if you change alignment settings, be sure you understand how unmapped/multi-mapped anchors are handled.

---

## 3. Legacy compatibility & engine design

Internally, `find_circ3` has **two engines**, selected automatically based on the **QNAME format** in the anchor SAM:

1. **Legacy engine (`legacy_call_iter`)**  
   - Activated when the SAM QNAMEs look like **unmapped2anchors3 output**, e.g.  
     `ERR2139486.12345_A__ACGT...` and `ERR2139486.12345_B`.
   - This path is a **close, line-by-line port of `find_circ.py`** to Python 3.
   - It uses the same breakpoint search, scoring, and classification logic.
   - It writes optional stats and spliced-read FASTA files, just like the original.

2. **Modern engine (HitAccumulator-based)**  
   - Activated for simple QNAMEs (e.g. `read1`, `read2`), such as the `tiny.sam` fixture.
   - Designed to be simpler and easier to extend while still reproducing the
     expected behavior for these tests.

The dispatcher lives in `run_find_circ`:

- If the SAM file contains `_A__`, `_A`, or `_B` style legacy anchor labels → **use legacy engine**.  
- Otherwise → **use modern engine**.

This means:

- For **real pipelines based on the original unmapped2anchors3**, you get
  behavior extremely close to legacy `find_circ.py`.
- For **test fixtures / plain SAM files**, you get behavior tuned for the
  `tests/test_against_legacy` expectations.

---

## 4. ERR2139486 chr21 regression: how close are we?

A key integration test is the chr21 subset:

```bash
TEST=tests/test_data/ERR2139486_chr21

time find-circ3 call \
  "$TEST/ERR2139486.anchors.sam" \
  --genome tests/data/ref/chr21.fa \
  --name ERR2139486_chr21 \
  --prefix ERR2139486_ \
  --anchor 20 \
  > "$TEST/ERR2139486_splice_sites.bed"
```

Compare against the original `find_circ.py` output using the helper script:

```bash
python2.7 tests/find_circ_py/cmp_bed.py \
  "$TEST/ERR2139486_splice_sites_legacy.bed" \
  "$TEST/ERR2139486_splice_sites.bed" \
  | tail
```

With the current implementation you should see:

- **53 / 54** junctions match exactly (coordinates + annotations).
- The **canonical circRNA** on chr21 is reproduced perfectly:
  - `chr21 33523435 33527305 ERR2139486_circ_000001 ... CIRCULAR,PERFECT_EXT,UNAMBIGUOUS_BP`
- Exactly **one** legacy junction is missing in the new output:

  ```text
  MISSING chr21 39127938 46691128 ERR2139486_circ_000002  25  ...
  ```

This is a **single ultra-long circular junction** spanning ~7.6 Mb.

### Why is that one circ missing?

Because the new CLI applies `--max-intron 200000` by default:

- span = `46691128 - 39127938 = 7,563,190 bp`  
- `7,563,190 > 200,000` → filtered as “implausibly large circ”.

If you want *strict, bit-for-bit legacy parity* for this dataset, you can simply relax the filter:

```bash
find-circ3 call \
  "$TEST/ERR2139486.anchors.sam" \
  --genome tests/data/ref/chr21.fa \
  --name ERR2139486_chr21 \
  --prefix ERR2139486_ \
  --anchor 20 \
  --max-intron 8000000 \
  > "$TEST/ERR2139486_splice_sites_full.bed"
```

Then `cmp_bed.py` will report **54 / 54** overlapping junctions.

In other words:

- **Default behavior** = legacy-like, but with a *sane intron-length guard* for real data.  
- **Parity mode** = raise `--max-intron` (or disable it) when you really want every long-range junction.

---

## 5. Recommended defaults & knobs

For most practical runs:

- `--anchor 20`  
  Match the anchor size used by your anchor-generation step.

- `--min-as-xs 2`  
  (Default) Use an AS–XS margin of 2 to mark anchors as “unique enough”.  
  This maps directly onto the legacy `min_uniq_qual` behavior.

- `--max-intron 200000`  
  Practical upper bound for intron length to avoid bizarre ultra-long candidates.  
  Increase only if you know what you’re doing.

- `--min-support 2` or higher  
  For noisy real-world data, requiring at least 2 supporting reads is often reasonable.

- `--allow-non-canonical`  
  Only enable if your biology / alignment strategy genuinely suggests non-GT/AG junctions are interesting; otherwise you’ll greatly increase ambiguity.

---

## 6. Running the test suite

Before and after major code changes, always run:

```bash
pytest -q
```

Key tests:

- `tests/test_against_legacy/test_regression.py`  
  - `test_cli_tiny_sam_matches_expected`  
  - `test_cli_cdr1as_emits_circular_junctions`

- `tests/test_breakpoints.py`  
  Unit tests for breakpoint logic and motif handling.

- `tests/test_anchors_cli.py`  
  Sanity checks for `find-circ3 anchors`.

- `tests/test_integration_bam_anchors_call.py`  
  An integration test currently marked **xfail** by design.

If these are green, you’re still on the rails and not falling into the rabbit hole.

---

## 7. Development notes / gotchas

1. **Legacy vs modern engine**
   - If you change how QNAMEs are parsed (`_A__`, `_A`, `_B`), re-run:
     - `tests/test_against_legacy/test_regression.py`
     - chr21 regression (`cmp_bed.py`) to ensure nothing broke.

2. **Read reconstruction**
   - The legacy engine tries to reconstruct the original read sequence in this order:
     1. QNAME-embedded sequence (`READ_A__SEQ`),
     2. `A.query_sequence`,
     3. `B.query_sequence`.
   - If you alter how anchors are named, make sure this logic still works.

3. **Anchors alignment**
   - Avoid casually adding `--no-unal` or other aggressive filters to the **anchor**-alignment step; it can silently drop informative anchors and make debugging harder.
   - Keep alignment parameters close to those used in tests when in doubt.

4. **New filters**
   - Any new filter (MAPQ thresholds, span limits, motif constraints, etc.) should be *off by default* or carefully documented and accompanied by regression updates.

5. **When in doubt: re-run chr21**
   - The chr21 subset (`ERR2139486_chr21`) plus `cmp_bed.py` is a practical regression test for real-world-like data.
   - If you refactor core logic and chr21 stays at “53/54 with only the 7.6 Mb circ missing by max-intron” (or “54/54” if you relax `--max-intron`), you’re still aligned with the original design.

---

## 8. Citation / provenance

This project reimplements and modernizes the original **find_circ.py**:

- Marvin Jens et al., 2013: the first **find_circ** circular RNA detection pipeline.

If you use `find_circ3` in a publication, please cite the original find_circ paper and, if appropriate, this repository.

---

Happy circ hunting, and may your anchors be unique and your junctions canonical (unless you really want them not to be). 🧬
