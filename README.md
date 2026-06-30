# RegCompassR

RegCompassR is an R package scaffold for RegCompass-Multiome: a multiome-supported, GPR-aware, sample-aware framework for reaction capacity potential and selected network-constrained feasibility analysis.

## v0.1-v0.5 scope

The v0.1 implementation focuses on the input layer for already annotated Seurat v4/Signac single-cell multiome objects. The v0.2 implementation adds sample-aware micropooling and pool-level pseudobulk summaries. The v0.3 implementation adds simple GPR parsing and Layer 1 reaction capacity potential. The v0.4 implementation adds pool/GPR/Q95 diagnostics with bootstrap uncertainty. The v0.5 implementation adds a minimal toy GEM plus baseline and selected demand QP wrappers around OSQP/rosqp. It does not rerun clustering, WNN, or flux/QP modeling.

Implemented functions:

- `rc_validate_seurat()` checks that a Seurat object contains the requested RNA assay, ATAC assay, required sample and cell-type metadata, optional condition/batch metadata, and optional embedding.
- `rc_extract_inputs()` validates the object and extracts RNA assay data, ATAC assay data, cell metadata, and an optional embedding into a plain R list.
- `rc_make_pools()` creates sample-aware micropools within sample, optional condition, cell type, and optional local-state/cluster strata without mixing cells across samples.
- `rc_pool_mean()` computes pool-level sparse means for normalized expression/residual matrices.
- `rc_pool_detection()` computes pool-level detection rates from raw counts for later dropout-aware correction.
- `rc_parse_gpr_simple()` / `rc_parse_gpr_table()` parse curated simple GPR rules.
- `rc_run_layer1_capacity()` computes GPR-aware Layer 1 reaction capacity and Q95 diagnostics.
- `rc_pool_diagnostics()` reports v0.4 pool-level diagnostics for depth, low-power pools, and metabolic/GPR gene detection.
- `rc_q95_bootstrap()` adds bootstrap confidence intervals for reaction-wise Q95 diagnostics.
- `rc_toy_gem()`, `rc_build_baseline_qp()`, `rc_solve_qp()`, and `rc_demand_qp()` provide the v0.5 toy GEM/QP MVP.

## Expected input

A Seurat v4 object should contain:

- RNA assay: commonly `RNA` or `SCT`
- ATAC/chromatin assay: commonly `ATAC` or `peaks`
- metadata columns such as `sample_id`, `cell_type`, optional `condition`, and optional `batch`
- optional reductions such as WNN, UMAP, PCA, LSI, or Harmony

## Example

```r
library(RegCompassR)

rc_validate_seurat(
  object,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition"
)

inputs <- rc_extract_inputs(
  object,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  celltype_col = "cell_type"
)
```

## v0.2 micropooling example

```r
pool_map <- rc_make_pools(
  inputs$meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition",
  state_col = "seurat_clusters",
  target_size = 80,
  min_size = 30,
  seed = 1
)
```

## v0.2 pseudobulk example

```r
rna_pool_mean <- rc_pool_mean(inputs$rna, pool_map)
rna_detection <- rc_pool_detection(inputs$rna, pool_map)
```

Use raw counts for detection rates. Use normalized data or residual matrices for expression scores; do not use imputed matrices in the main GPR formula.

## v0.3 Layer 1 GPR capacity example

```r
gpr_table <- data.frame(
  reaction_id = c("R_HEX1", "R_PFK", "R_LDH"),
  gpr = c("HK1 or HK2 or HK3", "PFKM and PFKL", "LDHA or LDHB")
)

layer1 <- rc_run_layer1_capacity(
  gpr_table = gpr_table,
  pool_expression = rna_pool_mean,
  pool_detection = rna_detection,
  promiscuity_mode = "sqrt",
  tau = 0.08
)

reaction_capacity_L1 <- layer1$reaction_capacity_L1
reaction_confidence <- layer1$reaction_confidence
q95_diagnostics <- layer1$q95_diagnostics
```

Layer 1 capacity is a reaction capacity potential, not a true flux estimate. The AND rule uses a Boltzmann-weighted average biased toward the minimum; it is not a LogSumExp soft minimum.

## v0.4 diagnostics example

```r
gpr_genes <- unique(unlist(layer1$parsed_gpr, use.names = FALSE))
pool_diag <- rc_pool_diagnostics(
  pool_map,
  rna_counts = inputs$rna,
  atac_counts = inputs$atac,
  state_col = "seurat_clusters",
  metabolic_genes = gpr_genes,
  gpr_genes = gpr_genes
)

q95_diag <- rc_q95_calibrate(
  layer1$reaction_capacity_raw,
  min_direct = 100,
  bootstrap = TRUE,
  B = 200
)$Q
```

## v0.5 toy GEM and QP example

```r
toy <- rc_toy_gem()
qp <- rc_build_baseline_qp(
  toy,
  penalty = rep(1, length(toy$reaction_id)),
  lambda = 1e-4,
  atpm_rxn = "ATPM",
  atpm_min = 1
)

# Requires rosqp.
base_solution <- rc_solve_qp(qp, settings = list(verbose = FALSE))
base_status <- rc_osqp_status(base_solution)

demand_qp <- rc_demand_qp(qp, reaction_id = "BIOMASS", delta = 1)
demand_solution <- rc_solve_qp(demand_qp, settings = list(verbose = FALSE))
```

v0.5 intentionally stays on the toy GEM / selected demand QP path and does not attempt full Human-GEM, all-reaction QP, FVA, thermodynamic constraints, or community exchange.
