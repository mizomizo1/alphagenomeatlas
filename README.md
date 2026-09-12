# AlphaGenome Atlas AVI scoring for Evo2 cohort

Reproducible code for scoring the Evo2 482-variant cohort with AlphaGenome Atlas AVI score and benchmarking it against Evo2 delta score.

This repository intentionally does **not** include:

- cohort data files
- AlphaGenome API keys
- generated output tables
- local notebooks

Keep those locally or in private storage.

## What this pipeline does

1. Reads a local `all_482_for_avi.tsv` file.
2. Uses the pre-harmonized `pos_hg38` coordinate for every variant.
3. Does **not** use the original hg19 coordinate for Atlas scoring.
4. Queries AlphaGenome Atlas `AVI_SCORE`.
5. Keeps both successful and failed variants in the output.
6. Summarizes AVI coverage by label and variant type.
7. Computes AUROC for AVI and Evo2 on available/matched subsets.
8. Runs a paired DeLong test for AVI vs Evo2 on the matched subset.
9. Runs leave-one-gene-out AUROC stability checks.

## Expected input columns

The scoring script expects a local TSV with at least these columns:

```text
record_id
label
set_type
gene
chrom
pos_hg19
pos_hg38
ref
alt
variant_type
evo2_delta_7b
spdi_hg38_0based
liftover_status
```

Notes:

- `label`: 1 = pathogenic, 0 = benign/control.
- `pos_hg38`: 1-based hg38 coordinate, used for AlphaGenome Atlas queries.
- `spdi_hg38_0based`: kept for audit. The main scoring script uses `chrom`, `pos_hg38`, `ref`, and `alt`.
- `ref` or `alt` may be `-` for indels; the main scoring script converts this to an empty allele and records any API failure.
- For indels, run the dedicated validation script below before concluding that the Atlas API/scorer lacks coverage.

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Set the API key as an environment variable. Do not put it in code or commit it.

```bash
export ALPHA_GENOME_API_KEY="your_api_key_here"
```

## Run scoring

```bash
python scripts/score_avi_all482.py \
  --input-tsv /path/to/all_482_for_avi.tsv \
  --outdir /path/to/alphagenome_avi_outputs
```

The script writes:

```text
01_all482_with_avi.tsv
02_failed_or_missing_avi.tsv
03_run_metadata.json
04_scorer_metadata.tsv
```

## Validate INDEL query representation

If INDELs are all missing, do not immediately conclude that AlphaGenome Atlas cannot score indels.
First test whether the query representation matches the Atlas representation.

```bash
python scripts/validate_indel_query_representations.py \
  --input-tsv /path/to/all_482_for_avi.tsv \
  --outdir /path/to/indel_validation_outputs \
  --use-ucsc-reference \
  --interval-probe
```

For stricter reproducibility, use a local hg38/GRCh38 FASTA instead of UCSC REST:

```bash
python scripts/validate_indel_query_representations.py \
  --input-tsv /path/to/all_482_for_avi.tsv \
  --outdir /path/to/indel_validation_outputs \
  --hg38-fasta /path/to/hg38.fa \
  --interval-probe
```

This tests, for each INDEL:

- direct `pos_hg38` + `ref/alt`
- direct `spdi_hg38_0based`
- off-by-one sanity checks
- left-padded previous-base representation
- `genome.normalize_variant(...)` representation when reference sequence is available
- small interval probes around the locus

It writes:

```text
indel_candidate_query_results.tsv
indel_candidate_strategy_summary.tsv
indel_interval_probe_results.tsv
indel_validation_run_metadata.json
```

## Run analysis

```bash
python scripts/analyze_avi_results.py \
  --scored-tsv /path/to/alphagenome_avi_outputs/01_all482_with_avi.tsv \
  --outdir /path/to/alphagenome_avi_outputs
```

The script writes:

```text
05_coverage_by_label_and_variant_type.tsv
06_label_coverage_test.tsv
07_predictor_auroc_self_coverage.tsv
08_delong_avi_vs_evo2.tsv
09_leave_one_gene_out_auc_stability.tsv
10_analysis_summary.json
```

## Recommended interpretation

This analysis should be interpreted as a full-cohort predictor comparison, not as a pathogenic-SNV-only correlation analysis.

Key checks before using in a manuscript:

- Were all variants queried using hg38 coordinates?
- Are failed variants retained in the output?
- Does AVI coverage differ by label?
- Is missingness explained by variant type, especially indels?
- Have indel failures been checked with normalized/left-padded representations?
- Is the AUROC comparison done on a matched subset when comparing AVI and Evo2?
- Are AlphaGenome package/API versions and scoring date recorded?

## Safety / reproducibility notes

- Never commit API keys.
- Never commit private cohort data.
- Never silently rescue reference mismatches by ref/alt swapping.
- Keep all failed rows for coverage-bias checks.
- Record package versions, scoring date, input file name, and query convention.
