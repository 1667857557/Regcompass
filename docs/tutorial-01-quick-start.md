# Tutorial 1: one-shot workflow

This is the shortest complete path from a paired-cell RNA+ATAC Seurat object to
condition-comparable reaction scores. Equations are in
[Tutorial 3](tutorial-03-mathematical-model.md). Public API:
[functions.md](functions.md).

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

## GEM and publication-bound medium

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
  scenario = "cantor2017_hplm",
  species = "human"
)
```

### Built-in scenario policy

The only current built-in scenario is:

| `scenario` | Published definition | Encoded scope |
|---|---|---|
| `"cantor2017_hplm"` | Cantor et al., *Cell* 2017; doi:10.1016/j.cell.2017.03.023 | HPLM components with exact published concentrations and direct one-to-one GEM exchange mapping. |

RegCompass deliberately does not retain the previous `physiologic`,
`normal_human_plasma`, `mouse_plasma`, RPMI/DMEM, single-nutrient challenge,
`minimal`, `compass_model_bounds`, or `permissive_all_exchange` identifiers as
built-in biological scenarios. Those labels either combined multiple evidence
sources, depended on manufacturer-only formulations, represented incomplete
implementations, or were technical model settings rather than published media.

For another published extracellular environment, supply exact reaction-level or
metabolite-level rows with `scenario = NULL`. Every row must contain
`reference_label` and a valid `reference_doi`:

```r
published_custom <- data.frame(
  medium_scenario_id = "published_environment_2024",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE,
  reference_label = "Author et al., Journal 2024",
  reference_doi = "10.xxxx/published.article",
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = NULL,
  species = "human",
  custom_medium = published_custom
)
```

Published concentrations are provenance; RegCompass does not convert them into
measured transporter fluxes. `uptake_scale` is fixed at `1`. Ambiguous
salt-to-free-ion conversions are omitted rather than approximated. The formal
contract and retained HPLM evidence are documented in
[Published medium scenarios](medium-presets.md).

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
    min_cells = 100L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      condition_mix = 0.5,
      condition_weight = "equal",
      outer_nfolds = 5L,
      inner_nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE
    )
  ),
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
  ),
  layer1_args = list(
    projection_component = "condition",
    comparison_support = "auto",
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
Stage 1 parallelism.

## Canonical interpretation

- `condition_full_oof` is the primary regulatory and metabolic penalty route.
- Jointly estimable edges form the common-support component selected by
  `comparison_support`.
- A non-estimable edge contributes a structural zero in that condition.
- A predictor equal to zero in every input cell remains represented without a
  fitted coefficient and contributes zero.
- Stage 2 builds one graph per cell type while all conditions of that cell type
  share the graph; condition is applied after graph clustering.
- `regulatory_alpha = 1` and `gpr_and_method = "min"` are canonical.
- One medium-specific structural model is reused across all conditions and
  metacells.

The workflow does not calculate depth matching, common-depth restriction, alpha
sensitivity, zero-support sensitivity, or link-saturation propagation.

## Inspect outputs

```r
result$grn$condition_fit_status
result$metacells$input_design
result$layer1$gene_projection_condition_full_oof
result$layer1$gene_projection_common_oof
result$microcompass$penalty_condition_full_oof
result$microcompass$penalty_common_oof
result$microcompass$penalty_condition_unique_increment
result$reaction_ranking
result$condition_contrast
```

Use [Tutorial 2](tutorial-02-stepwise-audit.md) for restartable stages,
[Tutorial 4](tutorial-04-targeted-reaction-remapping.md) for targeted reaction
extension, and [Tutorial 5](tutorial-05-condition-differential-analysis.md) for
condition statistics.