# Workflow

## Phase 1: strict-stratum processing

`rc_run_regcompass()` divides cells by:

```text
condition × sample × cell type
```

Each upstream worker performs one of two ATAC paths.

Without matching fragment files:

```text
SuperCell2 metacells with one fixed gamma
→ sum object ATAC peak raw counts by metacell membership
→ minimum-metacell filter
→ Pando GRN
```

With matching fragment files:

```text
SuperCell2 metacells with one fixed gamma
→ aggregate fragments by metacell membership
→ pseudobulk MACS2/MACS3 peak calling within the strict stratum
→ quantify newly called peaks with Signac::FeatureMatrix()
→ rebuild the metacell ChromatinAssay
→ minimum-metacell filter
→ TF-IDF and Pando GRN
```

Both paths then continue through:

```text
Pando-derived reaction confidence
→ reaction meta-module expansion
→ local FASTCORE completion
```

Different strata may produce different de novo peak sets. Their saved ATAC
matrices are imported by sparse peak-name union with zero fill. This union is
used for bookkeeping and cross-stratum matrices; Pando itself remains fitted
independently within each strict stratum.

Strata that produce fewer than the required metacells are marked
`skipped_too_few_metacells` and do not enter Pando, global calibration or
scoring. The global stage starts after all non-skipped retained strata succeed
and every biological sample still has at least one analyzable stratum. The
upstream worker pool is stopped before global processing.

## Phase 2: global calibration and shared GEM

All metacell GPR-gene logCPM matrices are merged. By default, each biological
sample receives equal total weight during gene scaling and reaction-wise Q95
calibration.

Optional `limma::removeBatchEffect()` correction occurs after logCPM merging and
before gene scoring. The preserved design should include biological variables
such as condition and cell type.

Locally completed reaction sets are unioned and deduplicated. The shared medium
is applied, and FASTCORE repairs only global core directions that remain
incomplete.

## Phase 3: directional scoring

All metacells use the same stoichiometric model. Each metacell has its own
penalty vector derived from expression capacity and Pando confidence.

A fresh Layer 2 worker pool evaluates:

```text
shared model × metacell
```

## Output order

```text
00_strata/
01_stratum_workflows/
02_global_layer1.rds
03_global_meta_modules.rds
regcompass_global_metacell_result.rds
```

For fragment-enabled strata, `01_metacells/` also contains aggregated fragment
files, fragment manifests, fragment-derived `atac_counts.rds`, and
`peaks/called_peaks.rds`.
