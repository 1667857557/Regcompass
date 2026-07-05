# RegCompassR

RegCompassR is a 1.0 framework for **multiome-supported, GPR-aware, sample-aware reaction capacity potential** analysis from annotated Seurat v4 RNA+ATAC objects, with a Layer 2 design for selected stoichiometrically constrained and balance-penalized flux-potential optimization.

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

RegCompassR 1.0 keeps true flux inference, enzyme activity inference, causal regulator discovery, full medium-constrained FBA, and thermodynamic flux estimation outside the package claims. Layer 2 uses Layer 1 capacity as priors, bounds modifiers, penalties, and reaction-selection evidence for constrained flux-potential optimization rather than as measured flux.

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

## RegCompassR 1.0 Layer 2 design

RegCompassR 1.0 extends the current Layer 1 workflow without reinterpreting Layer 1 as measured flux. The recommended positioning is:

> **multiome-supported, GPR-aware, sample-aware reaction capacity potential + selected stoichiometrically constrained flux-potential optimization**

```text
Layer 1: Seurat v4 RNA+ATAC -> reaction capacity potential
Layer 2: use Layer 1 output as GEM priors, bounds, penalties, and reaction-selection scores
Layer 2 output: stoichiometrically feasible / balance-penalized flux potential, not true flux
```

### What Layer 2 borrows from COMPASS and scFEA

Layer 2 should absorb COMPASS-style genome-scale metabolic modeling where useful: a full GEM structure, steady-state mass balance (`S v = 0`), LP optimization with reaction-specific objectives, and expression-supported penalty minimization. It should not run full COMPASS over every single cell and every Human-GEM reaction by default, nor should it treat expression-derived penalties as true enzyme activities or require measured media constraints.

Layer 2 should also absorb scFEA-style ideas: relaxed mass balance (`S v ≈ 0`), metabolite imbalance/stress outputs, expression priors constrained by network balance, and local/module-like subnetworks for speed. RegCompassR should use deterministic, interpretable LP/QP formulations rather than training a neural network in the default R package workflow.

### Layer 2 input contract

Layer 2 consumes the existing Layer 1 object fields:

```r
layer1$C_rel
layer1$C_raw
layer1$reaction_confidence
layer1$q95_diagnostics
layer1$gpr_diagnostics
layer1$pool_meta
layer1$parsed_gpr
```

`C_rel` remains a **relative reaction capacity potential**. It may be used as a reaction prior, a pool/profile-specific upper-bound modifier, an optimization penalty source, or a reaction-selection score; it is not a direct flux bound by itself.

### Recommended run units for the three core biological questions

| Biological question | Default run unit | Main Layer 2 result |
|---|---|---|
| Same broad cell class across conditions | `sample_celltype` | condition differential constrained capacity |
| Reaction status within one broad cell class | `sample_celltype` plus optional pool diagnostics | reaction status / feasibility / imbalance |
| One broad cell class versus other cell classes | `sample_celltype` | target-vs-other metabolic activation |

For condition comparisons, pools are within-sample denoising units, not independent biological replicates. The recommended workflow is Layer 1 on pools, aggregation to sample × cell type × condition profiles, Layer 2 optimization on those profiles, and sample-level condition statistics.

### GEM object

Layer 2 uses a unified GEM object:

```r
gem <- list(
  S = S,                         # sparse metabolite × reaction matrix
  reactions = reactions,         # reaction_id, name, lb, ub, reversible, subsystem
  metabolites = metabolites,     # metabolite_id, name, compartment
  gpr_table = gpr_table,
  exchange_rxns = exchange_rxns,
  demand_rxns = demand_rxns,
  model_source = "Human-GEM"
)
```

Validation should require sparse `S`, matching dimensions and names, unique reaction/metabolite IDs, `lb <= ub`, explicit exchange/demand annotations, and a diagnosable mapping between Layer 1 reaction IDs and GEM reaction IDs.

### Aggregating Layer 1 to profiles

For a sample × cell type × condition unit `u`, aggregate pools by default with medians:

```text
C_tilde[r,u]    = median(C_rel[r,p] for p in u)
Conf_tilde[r,u] = median(Conf[r,p] for p in u)
support_fraction[r,u] = supported finite pools / total pools
```

This preserves sample-level replication and avoids treating multiple pools from one sample as independent replicates.

### Mapping capacity to bounds, penalties, and priors

Layer 2 maps aggregated Layer 1 evidence to a bounded evidence score:

```text
q[r,u] = clamp(C_tilde[r,u] * Conf_tilde[r,u], q_min, 1)
```

Recommended defaults are `q_min = 0.05` and `alpha = 0.05`. Unsupported GPR reactions should be down-weighted and reported, not hard-deleted. Exchange, demand, sink, transport, and non-enzymatic reactions should not be controlled by GPR confidence in the same way as enzyme-catalyzed reactions.

Profile-specific bounds and penalties are:

```text
ub[r,u] = ub_base[r] * (alpha + (1 - alpha) * q[r,u])
lb[r,u] = -ub[r,u] for reversible reactions, otherwise 0
penalty[r,u] = -log(q[r,u] + epsilon)
v_prior[r,u] = v_scale * q[r,u]
```

### Core optimization modes

#### `hard_lp`

COMPASS/FBA-like strict steady-state feasibility:

```text
maximize c'v
subject to S v = 0 and lb_u <= v <= ub_u
```

The output is `stoich_feasible_capacity`. Infeasible reactions should be recorded in diagnostics rather than causing a hard pipeline failure.

#### `compass_penalty`

A two-step COMPASS-like LP. Step 1 maximizes a target reaction capacity. Step 2 fixes the target reaction to at least `omega * vmax` and minimizes expression-unsupported absolute flux using split variables. Recommended default: `omega = 0.95`.

Outputs include `reaction_penalty`, `compass_like_score`, and `target_flux_fraction`. A high penalty means the target is difficult to realize with expression-supported network routes; it should not be interpreted as high activity.

#### `relaxed_balance_qp` / `relaxed_balance_lp`

scFEA-like balance-penalized optimization with deterministic convex solvers. The QP form minimizes deviation from the Layer 1 prior plus a mass-balance penalty:

```text
minimize sum_i w[i,u] * (v_i - v_prior[i,u])^2 + lambda * ||S v||_2^2
subject to lb_u <= v <= ub_u
```

The LP form uses imbalance slack variables and absolute prior deviations. Outputs include `balanced_flux_potential`, `metabolite_imbalance`, `imbalance_score`, and `reaction_residual`. The default main Layer 2 result should be `sample_celltype × relaxed_balance_qp`, with `hard_lp` and `compass_penalty` as auxiliary diagnostics.

### Reaction selection and speed

RegCompassR 1.0 should avoid all profiles × all Human-GEM reactions × full FVA by default. The recommended approach is Layer 1 scoring, candidate reaction selection, stoichiometric-neighbor expansion, local sub-GEM construction, and LP/QP/FVA only on selected reactions.

```r
rc_select_reactions_layer1(
  layer1,
  mode = "union",
  top_n = 300,
  min_C_rel = 0.2,
  min_confidence = 0.3,
  include_neighbors = TRUE,
  neighbor_depth = 1,
  include_exchange = TRUE,
  include_demand = TRUE
)
```

Mandatory support reactions include selected reactions, one-step upstream/downstream reactions sharing metabolites, exchange/transport/demand/sink reactions, ATP/NAD/NADP/CoA cofactor-balance reactions, and user-specified pathway reactions. Recommended defaults are `top_n = 300`, `max_reactions = 1000`, `neighbor_depth = 1`, `q_min = 0.05`, `alpha = 0.05`, `omega = 0.95`, and `lambda_balance` sensitivity values of `1`, `10`, and `100`.

### Proposed Layer 2 modules and APIs

Planned files:

```text
R/gem_io.R
R/gem_validate.R
R/layer1_aggregate.R
R/reaction_selection.R
R/bounds.R
R/penalty.R
R/lp_solver.R
R/qp_solver.R
R/layer2_flux.R
R/fva_selected.R
R/layer2_stats.R
R/layer2_report.R
```

Planned exported functions:

```r
rc_read_humangem_model()
rc_read_sbml_model()
rc_validate_gem()
rc_match_layer1_to_gem()
rc_aggregate_layer1_profiles()
rc_select_reactions_layer1()
rc_make_subgem()
rc_capacity_to_bounds()
rc_capacity_to_penalty()
rc_build_layer2_problem()
rc_solve_hard_lp()
rc_solve_compass_penalty_lp()
rc_solve_relaxed_balance_qp()
rc_run_layer2_flux()
rc_selected_fva()
rc_compare_condition_within_celltype()
rc_compare_celltype_vs_others()
rc_summarize_reaction_status()
rc_write_layer2_report_md()
```

LP solvers should prefer optional `highs`, `glpk`, or `gurobi`; QP should use optional `osqp`. Gurobi must not be a hard dependency. All solver wrappers should return a common structure with `status`, `objective`, `flux`, `imbalance`, `runtime`, `solver`, and `diagnostics`.

### Recommended user-facing APIs

Condition comparison within one cell type:

```r
layer2_cond <- rc_run_layer2_flux(
  layer1 = layer1,
  gem = gem,
  run_unit = "sample_celltype",
  celltype_col = "cell_type",
  sample_col = "sample_id",
  condition_col = "condition",
  target_celltype = "Oligodendrocyte",
  mode = c("relaxed_balance_qp", "hard_lp", "compass_penalty"),
  selection = "layer1_union",
  max_reactions = 1000
)

res_cond <- rc_compare_condition_within_celltype(
  layer2 = layer2_cond,
  target_celltype = "Oligodendrocyte",
  condition_col = "condition",
  contrast = c("cKO", "Control")
)
```

Reaction status in one cell type:

```r
status_ol <- rc_summarize_reaction_status(
  layer2 = layer2_cond,
  celltype = "Oligodendrocyte",
  group_by = c("condition"),
  min_confidence = 0.3,
  max_imbalance = 0.2
)
```

Target cell type versus others:

```r
layer2_ct <- rc_run_layer2_flux(
  layer1 = layer1,
  gem = gem,
  run_unit = "sample_celltype",
  mode = "relaxed_balance_qp",
  selection = "layer1_union"
)

res_ct <- rc_compare_celltype_vs_others(
  layer2 = layer2_ct,
  target_celltype = "Oligodendrocyte",
  reference = "all_others",
  block_col = "sample_id",
  condition_col = "condition"
)
```

### Validation requirements

Mathematical tests should cover GEM validation, toy-network LP behavior (`A_ext -> A -> B -> C_ext`), COMPASS-like penalty behavior, relaxed balance behavior, and Layer 1 aggregation. Biological validation should report capacity change, stoichiometric feasibility change, imbalance change, and confidence change rather than a single score. External validation should prioritize matched metabolomics, stable isotope flux data, bulk metabolic enzyme proteomics, known pathway markers, and perturbation directions.

### Interpretation boundary

RegCompassR 1.0 may claim **multiome-supported, GPR-aware, sample-aware, stoichiometrically constrained reaction capacity potential**. It should not claim true single-cell flux, true enzyme activity, causal regulator inference, complete medium-constrained FBA, or thermodynamic flux estimation.

A conservative manuscript statement is:

> RegCompassR 1.0 uses multiome-derived GPR evidence to construct sample-aware reaction capacity priors, then performs selected stoichiometrically constrained and balance-penalized optimization to identify reaction programs that are both molecularly supported and network-feasible.
