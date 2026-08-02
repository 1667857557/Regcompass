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


Condition-aware Stage 1 requires **Pando >= 1.6.3**, native condition ABI 6,
and its high-dimensional memory contract. Budget-approved small targets retain
the dense fast path; high-dimensional targets use sparse residual validation
and an exact matrix-free Schur PCG refit without full predictor-square
allocation. RegCompass runs the native self-test on the master and up to two
workers before fitting and has no R fallback. An incompatible installation
stops immediately.

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

### Built-in biological scenarios

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

| Scenario | Authoritative composition source | Named override |
|---|---|---|
| `normal_human_plasma` | HPLM: *Cell* 2017 plus updated HPLM: *Cell Metabolism* 2021; Plasmax in *Science Advances* 2019 is validation only | none |
| `mouse_plasma` | absolute mouse plasma/interstitial-fluid metabolomics: *Nature* 2026; limited quantitative secondary values from Gardner and Stuart 2024 | none |
| `high_glucose` | identical HPLM 2017/2021 background | glucose 25 mM; Han 2015 |
| `low_glucose` | identical HPLM 2017/2021 background | glucose 1 mM; Han 2015 |
| `high_lactate` | identical HPLM 2017/2021 background | lactate 20 mM; San-Millan 2020 |
| `low_lactate` | identical HPLM 2017/2021 background | lactate 0.5 mM; Cho 2025 |
| `low_glutamine` | identical HPLM 2017/2021 background | glutamine 0.5 mM; Visagie 2015 Methods |

The five challenge scenarios use the same basal nutrient composition. Only the
named treatment row is changed, so a high-versus-low comparison does not also
compare unrelated RPMI, DMEM, or Plasmax backgrounds. Plasmax is retained as an
independent validation source and is not numerically averaged with HPLM.

Inspect composition and challenge provenance:

```r
unique(medium_scenarios[, intersect(c(
  "medium_scenario_id",
  "medium_background_id",
  "composition_primary_reference_doi",
  "composition_validation_reference_doi",
  "background_reference_doi",
  "background_validation_reference_doi",
  "challenge_reference_doi",
  "scenario_construction"
), colnames(medium_scenarios))])
```

These are modelling environments, not measured transporter fluxes.
Concentration-derived target caps remain explicit sensitivity assumptions and
are intersected with the original GEM directionality.

### Several built-in scenarios

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

### User-defined medium composition

Reaction-level bounds can be supplied directly:

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

Metabolite-level availability is also supported:

```r
custom_metabolites <- data.frame(
  metabolite_name = c("glucose", "glutamine", "lactate"),
  available = c(TRUE, TRUE, TRUE),
  concentration_mM = c(5, 0.55, 1.6),
  uptake_fraction = c(0.2, 0.275, 0.08),
  target_exchange_flag = c(TRUE, TRUE, TRUE),
  required_match = TRUE,
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = NULL,
  species = "human",
  custom_metabolites = custom_metabolites
)
```

Built-in and custom scenarios may be generated together by supplying a built-in
scenario vector plus `custom_medium` or `custom_metabolites`. Full references and
interpretation rules are in [Medium scenarios and evidence](medium-presets.md).

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
