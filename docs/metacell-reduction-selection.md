# Select the RNA and ATAC reductions used for metacells

RegCompass Stage 2 constructs multimodal SuperCell2 metacells from two cell-level embedding geometries:

- RNA reduction: `pca` dimensions `1:30` by default;
- ATAC reduction: `lsi` dimensions `2:30` by default;
- base random seed: `12345L` by default.

The RNA reduction is selectable. A precomputed Harmony reduction can replace PCA without changing how metacell RNA or ATAC counts are aggregated.

## Complete Stage 2 parameter example

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  sample_col = NULL,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30,
    seed = 12345L,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  )
)
```

The implementation defaults are:

```r
rna_reduction = "pca"
rna_dims = 1:30
atac_reduction = "lsi"
atac_dims = 2:30
gamma = 30L
seed = 12345L
min_cells_per_stratum = 100L
min_metacell_size = 20L
min_metacells_per_stratum = 2L
overwrite = FALSE
```

The same fields can be supplied to `rc_run_regcompass()` or `rc_run_regcompass_one_shot()` through `metacell_args`.

## Seed handling

`seed` is the reproducible base seed. RegCompass builds condition strata in a deterministic order and passes a stratum-specific seed to SuperCell2:

```text
seed_for_stratum = seed + stratum_index - 1
```

Thus repeated runs with the same cells, metadata, reductions, dimensions, and seed use the same seed sequence. Changing the ordering or identity of condition strata changes the assigned stratum-specific seeds and invalidates the cache contract.

## Use Harmony for the RNA geometry

```r
stopifnot(
  "harmony" %in% names(A@reductions),
  "lsi" %in% names(A@reductions),
  ncol(SeuratObject::Embeddings(A[["harmony"]])) >= 30,
  ncol(SeuratObject::Embeddings(A[["lsi"]])) >= 30
)

step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  sample_col = NULL,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "harmony",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30,
    seed = 12345L,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    overwrite = TRUE
  )
)
```

## What each geometry parameter changes

- `rna_reduction`: the RNA embedding used to calculate multimodal cell similarity.
- `rna_dims`: the exact RNA coordinates supplied to SuperCell2.
- `atac_reduction`: the ATAC embedding used to calculate multimodal cell similarity.
- `atac_dims`: the exact ATAC coordinates supplied to SuperCell2. LSI dimension 1 is excluded by default because it often tracks sequencing depth.
- `gamma`: approximate cells-per-metacell compression target.
- `seed`: deterministic base seed for SuperCell2.

The reduction names must exist in `A@reductions`, and the maximum requested dimension must not exceed the number of columns in the corresponding embedding.

## What Harmony changes

Harmony changes the RNA neighbourhood geometry used by SuperCell2 to decide which cells can be compressed together. It does not replace the RNA assay, alter RNA counts, or generate metacell expression values from Harmony coordinates. After memberships are determined, RegCompass aggregates the original RNA and ATAC assay counts.

This distinction matters biologically:

- use PCA when biological condition differences and technical batch differences are already adequately controlled;
- use Harmony when known technical batches distort the RNA geometry and the Harmony model was fitted without removing the biological contrast that RegCompass is intended to analyse;
- do not use a Harmony reduction that integrated away the condition effect of interest.

Condition remains the only hard metacell stratum. `celltype_col` remains the SuperCell2 construction label regardless of whether PCA or Harmony is selected.

## Cache and restart behaviour

The Stage 2 cache contract records:

- reduction names;
- selected dimensions;
- fingerprints of the selected embeddings;
- ordered cells and labels;
- RNA and ATAC assay fingerprints;
- `gamma`, `seed`, and metacell thresholds.

Changing any of these inputs invalidates existing Stage 2 checkpoints. Rebuild explicitly:

```r
metacell_args = list(
  rna_reduction = "harmony",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30,
  seed = 12345L,
  overwrite = TRUE
)
```

## Verify the selected reduction and seed

```r
step2$pooled$cache_contract$analysis_args[c(
  "rna_reduction",
  "rna_dims",
  "atac_reduction",
  "atac_dims",
  "gamma",
  "seed",
  "min_cells_per_stratum",
  "min_metacell_size",
  "min_metacells_per_stratum"
)]

step2$pooled$cache_contract$rna_reduction$embedding
step2$pooled$cache_contract$atac_reduction$embedding
```

The embedding objects are fingerprints of the exact coordinates used to construct the metacells.
