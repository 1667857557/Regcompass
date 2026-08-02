# Tutorial 1: one-shot workflow

This is the shortest complete path from a paired-cell RNA+ATAC Seurat object to
condition-comparable reaction scores. Public API: [functions.md](functions.md).
The condition-GRN equations and interface contract are described in
[condition-comparable-grn.md](condition-comparable-grn.md).

## Required object state

The object must contain paired RNA and ATAC assays for the same cells,
condition and broad-cell-type metadata, RNA PCA, ATAC LSI, and genome-compatible
peak coordinates.

```r
stopifnot(
  all(c("Group", "cell_type") %in% colnames(A@meta.data)),
  "pca" %in% names(A@reductions),
  "lsi" %in% names(A@reductions)
)
```

Stage 1 fixes the minimum retained broad-cell-type size at 300 paired cells. For
multiple conditions, each retained cell type must contain at least two eligible
condition strata. The condition model does not use sample/donor labels, nested
cross-validation, a sparse-group lambda path, or a native condition solver.

## GEM and medium

```r
library(RegCompassR)
library(BSgenome.Hsapiens.UCSC.hg38)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

Supported biological scenarios are:

```text
normal_human_plasma
mouse_plasma
high_glucose
low_glucose
high_lactate
low_lactate
low_glutamine
custom
```

Several built-in scenarios may be generated together:

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = c(
    "normal_human_plasma",
    "high_glucose",
    "low_glucose",
    "high_lactate",
    "low_lactate",
    "low_glutamine"
  ),
  species = "human"
)
```

Reaction-level custom bounds remain supported:

```r
custom_medium <- data.frame(
  medium_scenario_id = "my_measured_medium",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE,
  reference_label = "Optional experiment or publication label",
  reference_doi = "10.xxxx/optional.reference",
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

Full references and interpretation rules are in
[Medium scenarios and evidence](medium-presets.md).

## Run

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    k.knn = 30L,
    seed = 12345L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
  ),
  layer1_args = list(
    projection_component = "condition",
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
  upstream_workers = 6L,
  layer2_workers = 30L
)
```

Do not place `parallel` or `BPPARAM` inside `pando_infer_args`; the runner owns
Stage 1 parallelism. The following retired condition-GRN arguments are rejected:

```text
candidate_screen
condition_mix
condition_weight
alpha
nlambda / lambda / lambda_min_ratio
outer_nfolds / inner_nfolds
lambda_selection
scale
engine_control
comparison_conditions
```

## Condition-GRN algorithm

For every broad cell type with at least two eligible conditions:

1. Pando discovers candidate TF–peak–target edges on all eligible cells of that
   cell type.
2. Pando repeats candidate discovery separately in each condition.
3. Complete `(target, TF, region)` triples are unioned; TF, peak and target node
   sets are never recombined by Cartesian product.
4. The resulting target-specific dictionary is frozen.
5. Every condition fits the same Gaussian identity model
   `target ~ TF:peak`, using the same globally preprocessed RNA/ATAC layers and
   `scale = FALSE`.
6. Condition effects are the fitted coefficients. No global-coefficient
   calibration is applied.
7. BH-adjusted P values are computed across the fitted condition network. Only
   coefficients with `padj < 0.05`, sufficient absolute effect if requested, and
   target-model `R² >= min_model_rsq` enter the regulatory penalty.

The unfiltered coefficient, standard error, P value, adjusted P value,
estimability and rank diagnostics are retained in Stage 1 artifacts. An
unavailable coefficient remains `NA`; it is not silently interpreted as a fitted
zero. Ordinary GLM P values are conditional on the frozen candidate dictionary
and do not include selective-inference correction for candidate discovery.

## No condition or one condition

When `condition_col = NULL`, the column is absent, or only one non-missing level
is observed, Stage 1 directly runs the original Pando Gaussian interaction GRN
independently for each retained broad cell type. No condition coefficient or
condition fit contract is manufactured. The standard-Pando edge filter remains
`padj < 0.05` together with the configured model-fit and absolute-effect gates.

## Metacells and penalty handoff

Stage 2 builds one multimodal WNN graph per broad cell type. All conditions of
that type share the graph, and condition is applied after graph clustering so
final metacells remain condition-pure.

For multiple conditions, Pando reconstructs each retained paired-cell predictor
`TF RNA × peak ATAC`, multiplies it by the BH-significant `penalty_effect`, sums
by target, and only then aggregates by exact SuperCell membership. RegCompass
does not recompute TF×ATAC from metacell averages and does not refit or
renormalize coefficients after aggregation.

`regulatory_alpha = 1` and `gpr_and_method = "min"` remain canonical. The same
GEM, reaction order, bounds, direction and medium-specific `vmax` are reused
across conditions and metacells.

## Inspect outputs

```r
result$grn$condition_fit_status
result$grn$tf_peak_gene_condition_effect_all
result$grn$tf_peak_gene_condition_effect
result$metacells$input_design
result$layer1$gene_regulatory_modifier
result$microcompass$penalty
result$reaction_ranking
result$condition_contrast
```

For backward compatibility, several Stage 4/5 fields retain historical names
containing `_oof` or `common`. In condition mode their current meanings are:

```text
gene_projection_condition_full_oof = BH-filtered fixed-dictionary full-fit projection
gene_projection_common_oof         = compatibility alias of the same projection
gene_projection_condition_unique_oof = zero compatibility decomposition
penalty_condition_full_oof          = penalty derived from the primary projection
penalty_common_oof                  = compatibility alias
penalty_condition_unique_increment  = zero compatibility decomposition
```

These names do not imply that the current method performs OOF estimation or a
common-support decomposition.

Use [Tutorial 2](tutorial-02-stepwise-audit.md) for restartable stages,
[Tutorial 4](tutorial-04-targeted-reaction-remapping.md) for targeted reaction
extension, and [Tutorial 5](tutorial-05-condition-differential-analysis.md) for
condition statistics.
