# RegCompassR

RegCompassR is a simplified Layer 1 tool for **multiome-supported, GPR-aware, sample-aware reaction capacity potential** analysis from annotated Seurat v4 RNA+ATAC objects.

The package intentionally keeps the main analysis narrow:

```text
Seurat v4 RNA+ATAC counts
→ sample × cell_type × optional condition micropools
→ RNA raw-count pool sum
→ log2(CPM + 1)
→ robust z-score by gene across pools
→ sigmoid gene score
→ sqrt promiscuity correction
→ GPR capacity:
   AND = Boltzmann minimum-biased average, tau = 0.20
   OR  = sum across isoenzyme groups
→ cell-type Q95 continuous shrinkage
→ reaction confidence = median gene confidence × observed GPR gene fraction
→ C_raw, C_rel, reaction_confidence, minimal diagnostics
```

- `rc_validate_seurat()` checks that a Seurat object contains the requested RNA assay, ATAC assay, required sample and cell-type metadata, optional condition/batch metadata, and optional embedding.
- `rc_extract_inputs()` validates the object and extracts RNA assay data, ATAC assay data, cell metadata, and an optional embedding into a plain R list.
- `rc_make_pools()` creates sample-aware micropools within sample, optional condition, cell type, and optional local-state/cluster strata without mixing cells across samples.
- `rc_pseudobulk_counts()` sums raw counts by pool; `rc_filter_empty_pools()` removes zero-library pools; `rc_logcpm()` computes pool-level `log2(CPM + 1)` for Layer 1.
- `rc_pool_detection()` computes pool-level detection rates from raw counts for confidence/diagnostics only.
- `rc_parse_gpr_simple()` / `rc_parse_gpr_table()` parse curated simple GPR rules.
- `rc_run_layer1_capacity()` computes GPR-aware Layer 1 reaction capacity and Q95 diagnostics.
- `rc_pool_diagnostics()` reports v0.4 pool-level diagnostics for depth, low-power pools, and metabolic/GPR gene detection.
- `rc_q95_bootstrap()` adds bootstrap confidence intervals for reaction-wise Q95 diagnostics.
- `rc_toy_gem()`, `rc_build_baseline_qp()`, `rc_solve_qp()`, and `rc_demand_qp()` provide the v0.5 toy GEM/QP MVP.
- `rc_select_reactions()` and `rc_estimate_selected_demand_qp()` provide the v0.6 selected-demand QP planning layer.
- `rc_sample_aggregate()`, `rc_lm_by_reaction()`, and `rc_rank_regulators()` provide the v0.7 sample-level statistics and regulator candidate-prioritization layer.

## Expected input

An already annotated Seurat v4/Signac multiome object with:

- RNA assay counts, usually `RNA`
- ATAC/chromatin assay counts, usually `ATAC` or `peaks`
- identical RNA and ATAC cell barcodes
- metadata columns:
  - `sample_id`
  - `cell_type`
  - optional `condition`

## Main workflow

```r
library(RegCompassR)

inputs <- rc_extract_seurat_v4(
  object,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition"
)

pool_map <- rc_make_pools(
  inputs$meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition",
  target_size = 80,
  min_pool_size = 30,
  min_group_size = 30,
  seed = 1
)

layer1 <- rc_run_layer1_from_counts(
  gpr_table = gpr_table,
  rna_counts = inputs$rna_counts,
  pool_map = pool_map,
  stratum_col = "cell_type",
  tau = 0.20
)

RegCompassR Layer 1 uses the adjusted plan order: raw counts → pool sum → pool-level normalization. Do not average cell-level residuals, TF-IDF, or imputed expression for the main GPR capacity input.

```r
rna_pb <- rc_pseudobulk_counts(inputs$rna_counts, pool_map, fun = "sum")
pool_meta <- rc_build_pool_metadata(pool_map, inputs$meta)
filtered <- rc_filter_empty_pools(rna_pb, pool_meta)
rna_logcpm <- rc_logcpm(filtered$counts)
pool_meta <- filtered$pool_meta
rna_detection <- rc_pool_detection(inputs$rna_counts, pool_map)
rna_detection <- rna_detection[, colnames(rna_logcpm), drop = FALSE]
```

Detection rates are retained for confidence/diagnostics only and do not directly modify the gene score.

### Pool expression

Layer 1 uses only:

layer1 <- rc_run_layer1_capacity(
  gpr_table = gpr_table,
  pool_expression = rna_logcpm,
  pool_detection = rna_detection,
  pool_meta = pool_meta,
  stratum_col = "cell_type",
  promiscuity_mode = "sqrt",
  and_method = "boltzmann",
  tau = 0.20
)

where \(count_{g,p}\) is the raw-count sum of gene \(g\) in pool \(p\).

Cell-level Pearson residuals, TF-IDF, scaled data, and imputed matrices are not valid main inputs because their averages are not equivalent to pseudobulk count normalization.

### Gene score

```r
gpr_genes <- unique(unlist(layer1$parsed_gpr, use.names = FALSE))
pool_diag <- rc_pool_diagnostics(
  pool_map,
  rna_counts = inputs$rna_counts,
  atac_counts = inputs$atac_counts,
  state_col = "seurat_clusters",
  metabolic_genes = gpr_genes,
  gpr_genes = gpr_genes
)

q95_diag <- rc_q95_calibrate(
  layer1$reaction_capacity_raw,
  min_direct = 100,
  bootstrap = TRUE,
  B = 500
)$Q
```

with clipping to \([-6,6]\), then:

\[
s_{g,p}=\sigma(z_{g,p})
\]

### GPR capacity

Promiscuous genes are downweighted by:

\[
s'_{g,p}=s_{g,p}/\sqrt{N_{rxn}(g)}
\]

For AND groups:

\[
w_{g,p}=\frac{\exp(-s'_{g,p}/\tau)}{\sum_h \exp(-s'_{h,p}/\tau)}
\]

\[
C_k(r,p)=\sum_g w_{g,p}s'_{g,p}
\]

Default \(\tau=0.20\). Smaller values such as 0.08 behave closer to a hard minimum; larger values move toward a mean-like AND. RegCompassR reports a hard-min sensitivity flag but keeps one main result.

For OR groups:

\[
C_{raw}(r,p)=\sum_k C_k(r,p)
\]

### Q95 calibration

Within each cell type:

\[
Q_r=\rho_n Q_{r,celltype}+(1-\rho_n)Q_{r,global}
\]

\[
\rho_n=\frac{n}{n+80}
\]

\[
C_{rel}(r,p)=\min\left(1,\frac{C_{raw}(r,p)}{Q_r+\epsilon}\right)
\]

`C_rel` is a relative reaction capacity potential, not enzyme activity and not a flux bound.

### Reaction confidence

When only RNA detection is available, reaction confidence defaults to:

\[
Conf_{r,p}=median(Detection_{g,p}:g\in GPR_r)\times(1-missing\_gpr\_gene\_fraction)
\]

If a user supplies a gene-level confidence matrix, the same missing-GPR penalty is applied.

## Minimal diagnostics

The main workflow returns:

`rc_sample_aggregate()` aggregates pool-level reaction scores to biological sample × annotated cell-type medians, so differential analysis does not treat pools as independent biological replicates. `rc_lm_by_reaction()` fits simple reaction-wise sample-level linear models and reports BH-adjusted q-values within model terms. `rc_rank_regulators()` combines direct/adjusted association and support evidence into candidate regulator rankings only; rankings are not causal driver claims.

## Export and report helpers

```r
sample_matrix <- rc_sample_aggregate(layer1$reaction_capacity_L1, pool_meta)
rc_export_sample_matrix(sample_matrix, "output/sample_capacity.tsv")
rc_export_long_table(layer1$reaction_capacity_L1, "output/pool_capacity_long.tsv", value_col = "C_rel")
rc_write_report_md("output/layer1_report.md", q95_diagnostics = layer1$q95_diagnostics, gpr_diagnostics = layer1$gpr_diagnostics, confidence = layer1$reaction_confidence)
```
