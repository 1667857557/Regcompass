# RegCompassR

RegCompassR is a Layer 1 plus Layer 2 tool for **multiome-supported, GPR-aware, sample-aware reaction capacity potential** analysis from annotated Seurat v4 RNA+ATAC objects. RegCompassR 1.0 keeps Layer 1 as the evidence generator and adds a single Layer 2 main algorithm: selected-subnetwork, multiome/GPR-weighted, COMPASS-like two-step penalty LP.

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
→ reaction confidence:
   multiome/user gene confidence = GPR-aware AND/OR aggregation
   RNA-only detection = GPR-aware AND/OR aggregation
→ C_raw, C_rel, reaction_confidence, structured diagnostics
```

- `rc_validate_seurat()` / `rc_validate_seurat_v4()` check that a Seurat object contains the requested RNA assay, ATAC assay, required sample and cell-type metadata, optional condition/batch/state metadata, and optional embedding.
- `rc_extract_inputs()` / `rc_extract_seurat_v4()` validate the object and extract RNA assay data, ATAC assay data, cell metadata, and an optional embedding into a plain R list.
- `rc_check_metadata()` and `rc_write_input_summary()` create input QC summaries for sample/cell-type counts, missing metadata, optional condition-batch balance, and state provenance.
- `rc_drop_na_grouping()` removes cells with missing pooling labels; `rc_make_pools()` creates sample-aware micropools within sample, optional condition, cell type, and optional local-state/cluster strata without mixing cells across samples. If one major cell class is selected with no contrast, it runs as a single-class analysis by default; optionally, non-target cell types can be kept as an explicit control group. `rc_make_pool_seed_replicates()` repeats pooling over multiple seeds for sensitivity analysis.
- `rc_pseudobulk_counts()` sums raw counts by pool; `rc_filter_empty_pools()` removes zero-library pools; `rc_logcpm()` computes pool-level `log2(CPM + 1)` for Layer 1; `rc_check_pseudobulk_mapping()` spot-checks pool membership against pseudobulk columns.
- `rc_pool_detection()` computes pool-level detection rates from raw counts for confidence/diagnostics only; `rc_atac_pool_logcpm()` computes pooled ATAC accessibility logCPM for multiome confidence.
- `rc_parse_gpr_simple()` / `rc_parse_gpr_table()` parse curated simple GPR rules.
- `rc_run_layer1_capacity()` computes GPR-aware Layer 1 reaction capacity, Q95 calibration, reaction confidence, and diagnostics.
- `rc_reaction_confidence_gpr_aware()` aggregates gene confidence by GPR structure: AND = softmin/min/mean and OR = max/prob_or/sum_sqrtK.
- `rc_q95_calibrate()` and `rc_q95_bootstrap()` calibrate `C_raw` to `C_rel` and add Q95/bootstrap diagnostics.
- `rc_pool_diagnostics()` reports v0.4 pool-level diagnostics for depth, low-power pools, and metabolic/GPR gene detection.
- `rc_parallel_lapply()` and `rc_default_bpparam()` provide automatic BiocParallel-backed multi-core execution for expensive pool-, reaction-, model-, and bootstrap-level loops when BiocParallel is installed.
- `rc_sample_aggregate()`, `rc_sample_summary()`, `rc_export_sample_matrix()`, `rc_export_long_table()`, and `rc_write_report_md()` provide sample-aware summaries, table exports, and Markdown reporting.

RegCompassR 1.0 intentionally keeps standalone hard LP, scFEA-like relaxed balance QP/LP, selected FVA, thermodynamic constraints, and causal regulator discovery outside the runnable package scope. Hard LP is used only internally to compute Layer 2 `vmax` feasibility diagnostics.

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
  stratum_col = "cell_type",
  reaction_confidence_method = "gpr_aware",
  low_confidence_quantile = NULL
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
Q95 calibration, and reaction confidence. With multiome/user gene confidence,
reaction confidence is GPR-aware by default for both multiome gene confidence
and RNA-only detection. Legacy median mode is retained only for reproducibility.

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
  reaction_confidence_method = "gpr_aware",
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
  stratum_col = "cell_type",
  reaction_confidence_method = "gpr_aware",
  low_confidence_quantile = NULL
)
```

When `atac_counts` and `peak_gene_links` are supplied, the wrapper computes
RNA and ATAC pool percentiles within `stratum_col`, peak-gene link confidence,
discrete-null-corrected RNA/ATAC concordance, Fisher-shrunk positive RNA/ATAC
association, and nonnegative gene confidence before reaction confidence is
calculated. Single-pool strata have undefined percentiles (`NA`), not maximal
confidence.


### Single-cell-class and optional other-cell controls

If the analysis is intended for one major cell class and no contrast column is
provided, select that class with `target_celltype`. RegCompassR then pools only
that class and does not create an artificial contrast:

```r
pool_map <- rc_make_pools(
  inputs$meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  target_celltype = "T cell",
  target_size = 80,
  min_pool_size = 30,
  min_group_size = 30
)
```

If the desired design is target-vs-all-other-cells, opt in explicitly:

```r
pool_map <- rc_make_pools(
  inputs$meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  target_celltype = "T cell",
  include_other_celltypes_as_control = TRUE,
  target_contrast_label = "T cell",
  other_contrast_label = "other"
)
```

This second mode keeps the original labels in `original_cell_type`, recodes the
analysis `cell_type` to target vs `other`, and creates `celltype_contrast` when
no `condition_col` is supplied.

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
  reaction_confidence_method = "gpr_aware",
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

For OR groups, the default `or_method = "sum_sqrtK"` dampens isoenzyme-rich reactions:

\[
C_{raw}(r,p)=\frac{\sum_k C_k(r,p)}{\sqrt{K_r}}
\]

Set `or_method = "sum"` only when undampened isoenzyme summation is desired.

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

`reaction_confidence` is evidence for whether the RNA/multiome data support a
reaction's GPR rule. It is **not** flux, enzyme activity, or causal regulator
evidence.

Default multiome/user-confidence mode is GPR-aware:

- AND/enzyme-complex groups use `softmin` by default, so the weakest required
  subunit limits the group while remaining less brittle than a hard minimum.
- OR/isoenzyme groups use `max` by default, so one well-supported isoenzyme or
  complete complex can support the reaction.
- Missing required genes make that AND group incomplete. The reaction remains
  estimable if another OR alternative is complete; otherwise confidence is `NA`.
- Fixed `low_confidence_threshold = 0.25` is not used by this mode. Leave
  `low_confidence_quantile = NULL` to keep continuous scores, or set a quantile
  such as `0.10` to flag the lowest-scoring reactions within each confidence
  source.

```r
conf <- rc_reaction_confidence_gpr_aware(
  gpr_list = layer1$parsed_gpr,
  gene_confidence = gene_confidence,
  and_method = "softmin",
  or_method = "max",
  low_confidence_quantile = NULL
)
```

Key output columns include `reaction_confidence`, `n_and_groups_total`,
`n_and_groups_complete`, `complete_and_group_fraction`,
`best_and_group_observed_fraction`, `any_incomplete_gpr_group_flag`,
`reaction_unsupported_by_complete_gpr_flag`, and the deprecated alias
`missing_required_subunit_flag`.

The deprecated reproducibility mode `method = "legacy_median"` uses median aggregation:

\[
Conf_{r,p}=median(Evidence_{g,p}:g\in GPR_r)\times(1-missing\_gpr\_gene\_fraction)
\]

Select it explicitly only when reproducing older reports with `reaction_confidence_method = "legacy_median"`.

## Outputs and diagnostics

The main workflow returns `C_raw`, `C_rel`, `reaction_capacity_L1`,
`reaction_confidence`, `reaction_confidence_method`,
`reaction_confidence_source`, `reaction_confidence_summary`, `gene_confidence_components`, `q95_diagnostics`,
`gpr_diagnostics`, `tau_sensitivity`, `promiscuity_sensitivity`,
`and_method_sensitivity`, `capacity_long`, and `parsed_gpr`.
`rc_run_layer1_from_counts()` also returns `pool_meta`.

Q95 diagnostics include `q95_power_class`, bootstrap fields when requested, and
`all_missing_reaction_flag`. Reactions with all-`NA` `C_raw` remain `NA` in
`C_rel` and should be excluded from downstream rankings, along with reactions
that have `reaction_unsupported_by_complete_gpr_flag` or `q95_power_class == "very_low"`. Use `rc_filter_valid_reactions()` and `rc_rank_reactions()` for consistent downstream filtering/ranking.

`rc_sample_aggregate()` aggregates pool-level reaction scores to biological sample × annotated cell-type medians, so downstream analysis does not treat pools as independent biological replicates. `rc_sample_summary()` provides long-form median/IQR summaries by sample, cell type, and optional condition. Statistical modeling and regulator ranking are left to downstream project-specific analyses rather than treated as built-in causal inference.


## Function quick reference

| Function(s) | Purpose | Main notes |
|---|---|---|
| `rc_validate_seurat()`, `rc_validate_seurat_v4()` | Validate Seurat/Signac inputs | Checks assays, metadata, barcodes, and optional embedding. |
| `rc_extract_inputs()`, `rc_extract_seurat_v4()`, `rc_get_assay_counts()` | Extract plain count matrices and metadata | Use counts, not scaled/imputed values, for Layer 1. |
| `rc_check_metadata()`, `rc_write_input_summary()` | Input QC | Summarizes sample/cell-type/condition coverage and missing labels. |
| `rc_make_pools()`, `rc_make_pool_seed_replicates()`, `rc_filter_active_pool_map()` | Build sample-aware micropools | Supports single target-cell-class pooling by default and optional target-vs-other control pooling. |
| `rc_pseudobulk_counts()`, `rc_filter_empty_pools()`, `rc_logcpm()`, `rc_pool_mean()`, `rc_check_pseudobulk_mapping()` | Pool expression preprocessing | Sum raw counts, remove empty pools, then compute pool logCPM. |
| `rc_pool_detection()`, `rc_pool_diagnostics()` | Detection and pool QC | Detection supports confidence diagnostics; it does not alter capacity scores. |
| `rc_parse_gpr_simple()`, `rc_parse_gpr_table()`, `rc_metabolic_gpr_genes()` | GPR parsing and gene extraction | Supports simple AND/OR GPR rules and Human-GEM-style tables. |
| `rc_download_humangem_gpr_table()` | Download/convert Human-GEM GPRs | Optional convenience for generating a `reaction_id`/GPR table. |
| `rc_gene_score()`, `rc_gene_zscore()`, `rc_safe_scale()`, `rc_sigmoid()` | Gene scoring | Robust row scaling plus sigmoid transformation of pool logCPM. |
| `rc_and_capacity()`, `rc_or_capacity()`, `rc_reaction_capacity()`, `rc_reaction_capacity_one()` | GPR capacity calculation | AND defaults to Boltzmann bottleneck; OR capacity method is configurable. |
| `rc_q95_shrink()`, `rc_q95_calibrate()`, `rc_q95_bootstrap()` | Q95 calibration | Produces bounded `C_rel`; all-missing reactions stay `NA`. |
| `rc_percentile_by_stratum()`, `rc_concordance_null_correct()`, `rc_fisher_shrink()`, `rc_link_confidence()`, `rc_gene_confidence()` | Multiome gene confidence | Converts RNA/ATAC concordance and peak-gene links into gene-level confidence; `return_components = TRUE` exposes RA, detection, link, QC, and GPR-observed diagnostics. |
| `rc_recompute_signac_peak_gene_links()`, `rc_filter_peak_gene_links_to_gpr()`, `rc_atac_pool_logcpm()` | ATAC link utilities | Recompute/filter metabolic peak-gene links for multiome confidence. |
| `rc_reaction_confidence_gpr_aware()`, `rc_reaction_confidence()`, `rc_reaction_confidence_summary()` | Reaction confidence | GPR-aware aggregation is the default for multiome and RNA-only evidence; legacy median is reproducibility-only. |
| `rc_run_layer1_capacity()`, `rc_layer1_capacity()` | Core Layer 1 runner | Starts from pool expression and optional confidence inputs. |
| `rc_run_layer1_from_counts()`, `rc_run_layer1_from_seurat()` | End-to-end runners | Start from counts or Seurat object and return capacities, confidence, and diagnostics. |
| `rc_capacity_sensitivity()`, `rc_and_method_capacity_long()`, `rc_and_method_sensitivity()` | Sensitivity analysis | Compare tau, promiscuity, and AND aggregation choices. |
| `rc_filter_valid_reactions()`, `rc_rank_reactions()` | Filtering and ranking | Exclude all-missing, unsupported-GPR, very-low-Q95 reactions before ranking. |
| `rc_sample_aggregate()`, `rc_sample_summary()` | Sample-aware summaries | Aggregate pool-level results to biological sample/cell-type summaries. |
| `rc_export_sample_matrix()`, `rc_export_long_table()`, `rc_write_report_md()` | Export/report helpers | Write matrices, long tables, and Markdown summaries. |
| `rc_default_bpparam()`, `rc_parallel_lapply()` | Parallel execution | Uses BiocParallel when installed; otherwise deterministic sequential fallback. |

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

## RegCompassR 1.0 Layer 2: COMPASS-like two-step penalty LP

RegCompassR 1.0 uses one Layer 2 main algorithm: a **selected-subnetwork,
multiome/GPR-weighted, COMPASS-like two-step penalty LP**. The goal is not to
estimate true metabolic flux, enzyme activity, absolute flux, metabolite
production rate, or causal regulation. Instead, Layer 2 asks whether a target
reaction can be supported at low Layer1-derived multiome/GPR penalty while
respecting genome-scale stoichiometric steady-state constraints.

Layer 1 keeps its role as the source of reaction priors (`C_rel`),
GPR/multiome evidence (`reaction_confidence`), candidate reaction filtering, and
sample-aware aggregation. `C_rel` is still a relative reaction capacity
potential; it is not a flux bound and is not enzyme activity.

The primary Layer 2 outputs are:

```text
L2_compass_like_score
L2_compass_like_penalty
L2_vmax_internal
L2_feasible_flag
L2_solver_status
```

### Algorithm

For each sample × cell type unit by default, Layer 2 maps Layer 1 evidence to a
penalty:

```text
E[r,u] = max(C_rel[r,u], epsilon_C) * max(Conf[r,u], epsilon_Conf)
P[r,u] = min(-log(E[r,u] + epsilon), penalty_cap)
```

Defaults are `epsilon = 1e-6`, `epsilon_C = 1e-3`, `epsilon_Conf = 1e-3`, and
`penalty_cap = 20`. Exchange, demand, sink, and support reactions should use a
fixed low or medium-policy-specific penalty rather than a GPR penalty.

For each target reaction, Step 1 maximizes the target reaction under `S v = 0`
and GEM bounds to obtain `L2_vmax_internal`. This hard LP is retained only as an
internal feasibility diagnostic. If Step 1 is infeasible or `vmax <= 1e-8`, the
reaction is marked infeasible and receives score 0.

Step 2 requires `v_target >= omega * vmax` (default `omega = 0.95`) and minimizes
network-wide penalty-weighted absolute flux:

```text
minimize sum_i P[i,u] * |v_i|
subject to S v = 0, bounds, and v_target >= omega * vmax
```

The implementation uses positive/negative variable splitting to keep the
problem linear.

### Main function

```r
layer2 <- rc_run_layer2_compass_lp(
  layer1 = layer1,
  gem = gem,
  unit = "sample_celltype",
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition",
  selection_method = "auto",
  top_n = 300,
  min_C_rel = 0.15,
  min_confidence = 0.25,
  neighbor_depth = 1,
  max_subgem_reactions = 1000,
  omega = 0.95,
  solver = "highs",
  time_limit = 60
)
```

The returned object contains score, penalty, feasibility, solver status, internal
`vmax`, penalty components, selected reactions with selection reasons, sub-GEM
diagnostics, medium policy, unit metadata, and the method string
`"COMPASS-like two-step penalty LP"`.

### Reaction/sub-GEM selection

Layer 2 does not run all Human-GEM reactions by default. It selects candidate
reactions from high `C_rel`, high `reaction_confidence`, optional user-specified
reactions or pathways, one-hop shared-metabolite neighbors, and required
exchange/transport/demand/sink support reactions. Defaults are `top_n = 300`,
`min_C_rel = 0.15`, `min_confidence = 0.25`, `neighbor_depth = 1`, and
`max_subgem_reactions = 1000`. All-missing Layer 1 reactions, unsupported
complete GPR reactions, very-low Q95 reactions, and completely blocked
non-support reactions should be excluded from primary interpretation.

### Interpretation

Use wording such as:

- reaction capacity potential is higher;
- stoichiometrically plausible support is stronger;
- multiome/GPR evidence is more consistent with the metabolic network;
- COMPASS-like penalty is lower.

Avoid wording such as true flux is higher, enzyme activity is higher, metabolite
production rate is quantified, or causal regulation is proven.

Standalone hard LP, scFEA-like relaxed balance QP/LP, and selected FVA are not
implemented as RegCompassR 1.0 Layer 2 algorithms. Hard LP is only Step 1 of the
COMPASS-like two-step penalty LP, while relaxed balance and FVA are reserved for
possible future diagnostic extensions.
