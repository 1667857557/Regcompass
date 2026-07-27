# RegCompassR

RegCompassR 1.8.9 implements a shared-background regulatory–metabolic workflow for paired single-cell RNA+ATAC data.

## Canonical architecture

```text
all conditions within one cell type
→ one validated GREAT-domain Pando structural TF–peak–metabolic-gene universe
→ condition-aware TF×peak/target observability filter
→ one shared model edge universe across conditions
→ condition-balanced direct condition-theta elastic net
→ derived cross-condition GRN backbone + zero-sum deviations
→ condition-stratified full-size bootstrap reproducibility
→ condition-specific active sub-GRNs and metabolic targets
→ complete-GPR condition core reactions
→ ordered subsystem / KEGG–Reactome / master-Rhea expansion
→ merged biological reaction catalogue
→ one medium-specific union GEM reused by every condition and metacell
→ RNA+ATAC penalties
→ directional COMPASS-like LP scoring
→ compact final analysis tables
```

For edge \(e=(TF,peak,target)\), RegCompass directly estimates the
condition-specific coefficient \(\theta_{e,c}\). The reported backbone and
deviation are derived summaries:

\[
\beta_e=\frac1C\sum_c\theta_{e,c},
\qquad
\delta_{e,c}=\theta_{e,c}-\beta_e,
\qquad
\sum_c\delta_{e,c}=0.
\]

All conditions of one cell type use the same structural edge dictionary,
filtered model edge universe, predictor scale, penalty structure, and candidate
ordering. The elastic-net L1 penalty acts directly on \(\theta_{e,c}\), so an
edge may be exactly zero in one condition and non-zero in another. A reaction
becomes a condition core only when at least one complete GPR branch is contained
in the condition target-gene set.

## Installation

The validated default profile remains Seurat v4:

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github(
  "1667857557/SuperCell_Seurat_V4@c8b94949cd8a5ff7403f9f186c516f8efbac9b6f"
)
remotes::install_github(
  "1667857557/Pando_regcompass@6f42c8143bec6610b001e714a51627337f6d9ba9"
)
remotes::install_github("1667857557/Regcompass")
```

Pando 1.1.3 or later is required. RegCompass validates the Pando version-2 design fingerprint before fitting the multitask GRN.

SeuratObject/Seurat 5.x with Signac 1.12–1.x is also accepted. See [Seurat compatibility](docs/seurat-compatibility.md).

## Required metadata and current SuperCell2 contract

The canonical workflow requires only:

```text
condition_col
celltype_col
```

No biological-sample column is accepted or interpreted by `rc_run_regcompass()`, `rc_regcompass_step_grn()`, `rc_regcompass_step_metacells()`, or `rc_make_supercell2_metacells()`.

Stage 2 follows the current SuperCell2 API:

```text
RegCompass splits the Seurat object by `strata_cols = condition_col`
→ calls SCimplify_for_Seurat(label = celltype_col)
→ validates exact label-pure membership
```

There is no artificial condition-pool metadata field and no `sample_col` adapter.

## Minimal complete run

```r
library(RegCompassR)
library(BSgenome.Hsapiens.UCSC.hg38)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem,
  scenario = "high_glucose",
  species = "human"
)

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "Group",
  celltype_col = "cell_type",

  grn_mode = "multitask_shared_backbone",
  pando_args = list(
    min_cells = 100,
    pando_design_args = list(
      peak_to_gene_method = "GREAT",
      upstream = 100000,
      downstream = 0,
      extend = 1000000,
      only_tss = FALSE,
      min_tf_detection = 0,
      min_peak_detection = 0,
      min_target_detection = 0,
      max_edges_per_target = Inf
    )
  ),
  multitask_args = list(
    alpha = 0.5,
    global_penalty_factor = 1,
    deviation_penalty_factor = 1,
    lambda_rule = "lambda.1se",
    nfolds = 5,
    n_bootstrap = 100,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    min_bootstrap_success_fraction = 0.8,
    min_cv_rsq = 0,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    min_detected_cells_per_condition = 10,
    min_detection_fraction_per_condition = 0.01,
    seed = 12345L
  ),

  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30,
    seed = 12345L,
    min_cells_per_stratum = 300,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  ),

  layer1_args = list(
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),

  medium_scenarios = medium_scenarios,
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      completion_time_limit = 600,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000,
      strict = TRUE
    )
  ),
  upstream_workers = 6,
  layer2_workers = 30
)
```

The canonical peak-to-gene rule is `GREAT`. With `extend = 1000000`, it creates
a broad distal regulatory-domain hypothesis without using target-expression
correlation to admit candidates. The subsequent condition-aware observability
filter and regularised direct-theta model determine which candidates have usable
support. `Signac` remains available only as an explicit sensitivity override.

The canonical default is `n_bootstrap = 100L`. At a true selection probability of 0.7, the binomial Monte Carlo standard error is approximately 0.046; RegCompass also reports Wilson 95% intervals for each edge.

## Grounded GRN parameter policy

The Pando structural design remains condition agnostic. Its pooled detection thresholds stay at zero so a regulator that is restricted to one condition is not removed before the multitask model. RegCompass then applies a condition-aware observability rule. For condition \(c\),

\[
m_c=\min\left(n_c,\max\left(10,\left\lceil0.01n_c\right\rceil\right)\right).
\]

An edge enters the shared model universe when, in at least one condition,

\[
\sum_{u\in c}I(T_{t,u}>0\land A_{p,u}>0)\ge m_c
\]

and

\[
\sum_{u\in c}I(Y_{g,u}>0)\ge m_c.
\]

The TF and peak must therefore produce a non-zero interaction predictor in the same cells. This filter uses detection only, not target correlation or fitted effect size.

`alpha = 0.5` retains both lasso sparsity and ridge stabilization for correlated TF–peak predictors. `global_penalty_factor` and `deviation_penalty_factor` are compatibility aliases for one common direct-theta penalty and must be equal. They no longer represent separately penalised backbone and deviation coordinates.

`min_sign_stability = 0.8` means at least 90% agreement on one sign among bootstrap fits in which the edge is selected, because \(\rho=|2q-1|\). An active edge must also have strictly positive out-of-fold \(R^2\), so it improves on the condition-centred null predictor.

Finite `max_edges_per_target` values are rejected in the canonical model. Pando candidate order is deterministic but is not an evidence ranking; arbitrary top-K truncation would change the biological hypothesis without a statistical criterion.

See [GRN parameter policy](docs/grn-parameter-policy.md) and [Shared-background multitask GRN mathematics](docs/multitask-shared-grn.md).

## Core calculations

The Pando-style candidate predictor is

\[
x_{e,u}=T_{t,u}A_{p,u}.
\]

Target expression and predictors are centred within condition. Edge scales are shared across conditions, and observation weights make each condition contribute the same total regression loss. The direct condition-theta objective is

\[
\min_\Theta
\frac12\sum_u w_u
\left(y_u^\circ-\sum_e\widetilde x_{e,u}\theta_{e,c(u)}\right)^2
+\lambda\alpha p_\theta\sum_{e,c}|\theta_{e,c}|
+\frac{\lambda(1-\alpha)p_\theta}{2}\sum_{e,c}\theta_{e,c}^2.
\]

For successful bootstrap fits:

\[
\Pi_{e,c}=
\frac{1}{B_s}\sum_{b\in\mathcal B_s}
I\!\left(|\theta_{e,c}^{(b)}|>\varepsilon\right),
\]

\[
\rho_{e,c}=
\left|
\frac{\sum_b I_e^{(b)}\operatorname{sign}(\theta_{e,c}^{(b)})}
     {\sum_b I_e^{(b)}}
\right|.
\]

The Layer 1 projection weight is

\[
\widetilde\theta_{e,c}
=\widehat\theta_{e,c}\Pi_{e,c}\rho_{e,c}.
\]

The full-size bootstrap estimates cell-resampling reproducibility. It is not described as formal half-sample stability-selection PFER control.

Layer 1 collapses TF coefficients sharing the same measured peak before ATAC projection. One target-specific normalization is shared across conditions. If no active edge exists, the modifier is zero and multiome support equals RNA-only support.

## Compact final result

Stage 6 deliberately does not embed full Stage 1–4 objects. Primary tables are:

```r
result$table_manifest
result$reaction_ranking
result$condition_contrast
result$active_regulatory_edges
result$condition_target_genes
result$core_reactions
result$meta_module_summary
result$grn_metacell_group_coverage
result$reaction_catalog
result$reaction_evidence
```

`reaction_ranking` and `condition_contrast` contain only comparison fields. Reaction name, formula, GPR, and cross-references are stored once per reaction in `reaction_catalog`.

`result$microcompass` remains available because downstream reaction tests, direction reports, and plots require unit-level Layer 2 scores. Detailed candidate edges, all coefficients, metacell matrices, Layer 1 matrices, and full module membership remain in their stage checkpoints:

```r
result$stage_provenance$detailed_sources
```

This reduces result size without arbitrary top-N filtering or discarding scored reactions.

## Inspectable stages

- `rc_regcompass_step_grn()`: structural candidates, observability-filtered model universe, multitask coefficients, CV, and bootstrap.
- `rc_regcompass_step_metacells()`: condition-stratified, cell-type-labelled SuperCell2 metacells.
- `rc_regcompass_step_meta_modules()`: target genes, complete-GPR cores, and biological modules.
- `rc_regcompass_step_layer1()`: RNA support, ATAC modifier, and reaction capacity.
- `rc_regcompass_step_layer2()`: shared union GEM and directional LP scoring.
- `rc_regcompass_step_results()`: compact analysis tables and provenance.
- `rc_regcompass_step_target_union()`: direct-equivalent second-pass targets in the cached model.

## Legacy mode

Independent condition-by-cell-type Pando fitting remains available only for reproducibility:

```r
result_legacy <- rc_run_regcompass_one_shot(
  ...,
  grn_mode = "legacy_condition_pando",
  pando_args = list(
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    ),
    padj_threshold = 0.05,
    min_model_rsq = 0.1
  )
)
```

Legacy adjusted-P-value thresholds are not interpreted as significance thresholds in multitask mode.

## Tutorials

1. [Quick start](docs/tutorial-01-quick-start.md)
2. [Stepwise audit](docs/tutorial-02-stepwise-audit.md)
3. [Restart and sensitivity](docs/tutorial-03-advanced-restart.md)
4. [Targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
5. [Condition analysis from GRN to reaction](docs/tutorial-05-condition-differential-analysis.md)

Additional references:

- [GRN parameter policy](docs/grn-parameter-policy.md)
- [Stage contracts](docs/stage-interface-contracts.md)
- [Medium presets](docs/medium-presets.md)
- [Direction-aware reporting](docs/direction-aware-condition-reporting.md)
