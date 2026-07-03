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

- `rc_validate_seurat()` / `rc_validate_seurat_v4()` check that a Seurat object contains the requested RNA assay, ATAC assay, required sample and cell-type metadata, optional condition/batch/state metadata, and optional embedding.
- `rc_extract_inputs()` / `rc_extract_seurat_v4()` validate the object and extract RNA assay data, ATAC assay data, cell metadata, and an optional embedding into a plain R list.
- `rc_check_metadata()` and `rc_write_input_summary()` create input QC summaries for sample/cell-type counts, missing metadata, optional condition-batch balance, and state provenance.
- `rc_drop_na_grouping()` removes cells with missing pooling labels; `rc_make_pools()` creates sample-aware micropools within sample, optional condition, cell type, and optional local-state/cluster strata without mixing cells across samples; `rc_make_pool_seed_replicates()` repeats pooling over multiple seeds for sensitivity analysis.
- `rc_pseudobulk_counts()` sums raw counts by pool; `rc_filter_empty_pools()` removes zero-library pools; `rc_logcpm()` computes pool-level `log2(CPM + 1)` for Layer 1; `rc_check_pseudobulk_mapping()` spot-checks pool membership against pseudobulk columns.
- `rc_pool_detection()` computes pool-level detection rates from raw counts for confidence/diagnostics only; `rc_atac_pool_logcpm()` computes pooled ATAC accessibility logCPM for multiome confidence.
- `rc_parse_gpr_simple()` / `rc_parse_gpr_table()` parse curated simple GPR rules.
- `rc_run_layer1_capacity()` computes GPR-aware Layer 1 reaction capacity and Q95 diagnostics.
- `rc_pool_diagnostics()` reports v0.4 pool-level diagnostics for depth, low-power pools, and metabolic/GPR gene detection.
- `rc_q95_bootstrap()` adds bootstrap confidence intervals for reaction-wise Q95 diagnostics.
- `rc_parallel_lapply()` and `rc_default_bpparam()` provide automatic BiocParallel-backed multi-core execution for expensive pool-, reaction-, model-, and bootstrap-level loops when BiocParallel is installed.
- `rc_sample_aggregate()`, `rc_sample_summary()`, `rc_export_sample_matrix()`, `rc_export_long_table()`, and `rc_write_report_md()` provide sample-aware summaries, table exports, and Markdown reporting.

RegCompassR v0.9 intentionally keeps GEM/QP solving, selected-demand QP planning, FVA, thermodynamic constraints, and causal regulator discovery outside the runnable package scope.

## Expected input

An already annotated Seurat v4/Signac multiome object with:

- RNA assay counts, usually `RNA`
- ATAC/chromatin assay counts, usually `ATAC` or `peaks`
- identical RNA and ATAC cell barcodes
- metadata columns:
  - `sample_id`
  - `cell_type`
  - optional `condition`


## Human-GEM GPR tables and metabolic peak-gene links

RegCompassR can download the official Human-GEM repository archive and convert its
model GPR rules into the simple `reaction_id`, `and_group_id`, `gene` table used
by the Layer 1 workflow:

```r
hg <- rc_download_humangem_gpr_table(
  destdir = "data/Human-GEM",
  ref = "main",
  gene_format = "symbol"
)

gpr_table <- hg$gpr_table
metabolic_genes <- hg$metabolic_genes
```

The returned `metabolic_genes` vector is the target gene set for recomputing
metabolic peak-gene links. RegCompassR now provides
`rc_recompute_signac_peak_gene_links()`, which calls `Signac::LinkPeaks()`
internally with `genes.use = metabolic_genes`, extracts `Signac::Links()` from
the ATAC assay, converts the result to the `peak_id`, `gene`, `weight` table,
and filters it to GPR metabolic genes. Signac remains an optional dependency, so
install Signac before running this step.

```r
# Preferred Seurat entry point: recompute metabolic Signac links by default.
layer1 <- rc_run_layer1_from_seurat(
  gpr_table = gpr_table,
  object = object,
  pool_map = pool_map,
  pool_meta = pool_meta,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  recompute_peak_gene_links = TRUE
)

# Lower-level equivalent when you want the link table explicitly.
peak_gene_links <- rc_recompute_signac_peak_gene_links(
  object = object,
  metabolic_genes = metabolic_genes,
  peak_assay = "ATAC",
  expression_assay = "RNA"
)

# The Signac-updated object is retained for users who want to inspect links in Seurat.
object <- attr(peak_gene_links, "seurat_object")

layer1 <- rc_run_layer1_from_counts(
  gpr_table = gpr_table,
  rna_counts = inputs$rna_counts,
  pool_map = pool_map,
  pool_meta = pool_meta,
  atac_counts = inputs$atac_counts,
  peak_gene_links = peak_gene_links,
  stratum_col = "cell_type"
)
```


## Parallel execution

RegCompassR parallelizes the most expensive embarrassingly parallel work units: pool pseudobulk/detection summaries, reaction-wise GPR capacity calculations, Q95 bootstrap diagnostics, pool diagnostics, and reaction-wise linear models. These functions accept `BPPARAM` for an explicit `BiocParallelParam`; when `BPPARAM = NULL`, they call `rc_default_bpparam()` and use all detected CPU cores except one when BiocParallel is installed. If BiocParallel is unavailable, or if only one worker is requested, the same code falls back to deterministic sequential `lapply()` execution.

Control the automatic worker count globally with either:

```r
options(RegCompassR.workers = 4)
# or before starting R:
# Sys.setenv(REGCOMPASS_WORKERS = "4")
```

Set `options(RegCompassR.workers = 1)` or pass `BPPARAM = FALSE` to force sequential execution for debugging or constrained environments.

## Main workflow

RegCompassR keeps the capacity input as raw RNA counts until after pooling. The
recommended wrapper is `rc_run_layer1_from_counts()`, which performs the same
calculation as the manual steps below: pseudobulk pool sums, pool-level
`log2(CPM + 1)`, robust gene z-scores, sigmoid gene scores, GPR capacity,
Q95 calibration, and reaction confidence.

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

pool_meta <- rc_build_pool_metadata(pool_map, inputs$meta)

layer1 <- rc_run_layer1_from_counts(
  gpr_table = gpr_table,
  rna_counts = inputs$rna_counts,
  pool_map = pool_map,
  pool_meta = pool_meta,
  stratum_col = "cell_type",
  promiscuity_mode = "sqrt",
  and_method = "boltzmann",
  tau = 0.20,
  bootstrap = TRUE,
  B = 500
)
```

To enable multiome-supported gene confidence, pass ATAC counts and curated peak-gene links:

```r
layer1 <- rc_run_layer1_from_counts(
  gpr_table = gpr_table,
  rna_counts = inputs$rna_counts,
  pool_map = pool_map,
  pool_meta = pool_meta,
  atac_counts = inputs$atac_counts,
  peak_gene_links = peak_gene_links,
  stratum_col = "cell_type"
)
```

When `atac_counts` and `peak_gene_links` are supplied, the wrapper computes
RNA and ATAC pool percentiles within `stratum_col`, peak-gene link confidence,
discrete-null-corrected RNA/ATAC concordance, Fisher-shrunk positive RNA/ATAC
association, and nonnegative gene confidence before reaction confidence is
calculated. Single-pool strata have undefined percentiles (`NA`), not maximal
confidence.

### Manual calculation flow

These steps are useful when inspecting intermediate objects or providing custom
gene confidence.

```r
rna_pb <- rc_pseudobulk_counts(inputs$rna_counts, pool_map, fun = "sum")
filtered <- rc_filter_empty_pools(rna_pb, pool_meta)
rna_logcpm <- rc_logcpm(filtered$counts)
pool_meta <- filtered$pool_meta

rna_detection <- rc_pool_detection(inputs$rna_counts, pool_map)
rna_detection <- rna_detection[, colnames(rna_logcpm), drop = FALSE]

layer1 <- rc_run_layer1_capacity(
  gpr_table = gpr_table,
  pool_expression = rna_logcpm,
  pool_detection = rna_detection,
  pool_meta = pool_meta,
  stratum_col = "cell_type",
  promiscuity_mode = "sqrt",
  and_method = "boltzmann",
  tau = 0.20,
  bootstrap = TRUE,
  B = 500
)
```

Detection rates are retained for confidence/diagnostics only and do not directly modify the gene score. Cell-level Pearson residuals, TF-IDF, scaled data, and imputed matrices are not valid main inputs because their averages are not equivalent to pseudobulk count normalization.

### Gene score

For pool-level logCPM expression, `rc_gene_score()` computes robust row-wise
z-scores clipped to `[-6, 6]`, then applies the logistic sigmoid:

\[
s_{g,p}=\sigma(z_{g,p})
\]

where `count_{g,p}` is the raw-count sum of gene `g` in pool `p` before logCPM normalization.

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

Within each `stratum_col` group:

\[
Q_r=\rho_n Q_{r,stratum}+(1-\rho_n)Q_{r,global}
\]

\[
\rho_n=\frac{n}{n+80}
\]

\[
C_{rel}(r,p)=\min\left(1,\frac{C_{raw}(r,p)}{Q_r+\epsilon}\right)
\]

`rc_q95_calibrate()` reports an ordered `q95_power_class` (`very_low`, `low`, `moderate`, `adequate`, `high`) and flags unstable Q95 estimates when `q95_ci_width / max(q_value, 1e-6) > 0.5`.

`C_rel` is a relative reaction capacity potential, not enzyme activity and not a flux bound.

### Reaction confidence

When only RNA detection is available, reaction confidence defaults to:

\[
Conf_{r,p}=median(Detection_{g,p}:g\in GPR_r)\times(1-missing\_gpr\_gene\_fraction)
\]

If multiome confidence or a user-supplied gene-level confidence matrix is available, the same missing-GPR penalty is applied to the median gene confidence.

## Minimal diagnostics

The main workflow returns `C_raw`, `C_rel`, `reaction_capacity_L1`,
`reaction_confidence`, `q95_diagnostics`, `gpr_diagnostics`,
`tau_sensitivity`, `promiscuity_sensitivity`, `and_method_sensitivity`,
`capacity_long`, and `parsed_gpr`. The `rc_run_layer1_from_counts()` wrapper also
returns `pool_meta` and `reaction_confidence_source`.

`rc_sample_aggregate()` aggregates pool-level reaction scores to biological sample × annotated cell-type medians, so downstream analysis does not treat pools as independent biological replicates. `rc_sample_summary()` provides long-form median/IQR summaries by sample, cell type, and optional condition. Statistical modeling and regulator ranking are left to downstream project-specific analyses rather than treated as built-in causal inference.

## Export and report helpers

```r
input_summary <- rc_check_metadata(inputs$meta, condition_col = "condition")
rc_write_input_summary(input_summary, "output/input_qc")

sample_matrix <- rc_sample_aggregate(layer1$C_rel, layer1$pool_meta)
sample_summary <- rc_sample_summary(layer1$C_rel, layer1$pool_meta, condition_col = "condition")
rc_export_sample_matrix(sample_matrix, "output/sample_capacity.tsv")
rc_export_long_table(layer1$C_rel, "output/pool_capacity_long.tsv", value_col = "C_rel")
rc_write_report_md(
  "output/layer1_report.md",
  q95_diagnostics = layer1$q95_diagnostics,
  gpr_diagnostics = layer1$gpr_diagnostics,
  confidence = layer1$reaction_confidence
)
```
