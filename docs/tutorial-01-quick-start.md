# Tutorial 1: one-shot workflow

This is the complete one-shot route from a paired-cell RNA+ATAC Seurat object
to condition-comparable reaction scores. The condition-GRN equations are in
[Tutorial 3](tutorial-03-mathematical-model.md), while restartable execution is
covered in [Tutorial 2](tutorial-02-stepwise-audit.md).

## Required object state

The input must contain paired RNA and ATAC assays for the same cells, broad
cell-type metadata, RNA PCA, ATAC LSI, and genome-compatible peak coordinates.
A condition column is optional: two or more effective levels activate the
common-dictionary condition-GRN route; no condition or one effective level uses
direct per-cell-type Pando GRNs.

```r
stopifnot(
  "cell_type" %in% colnames(A@meta.data),
  "pca" %in% names(A@reductions),
  "lsi" %in% names(A@reductions)
)
```

Stage 1 fixes the minimum retained broad-cell-type size at 300 paired cells. In
multi-condition mode, each fitted cell type must retain at least two eligible
condition strata.

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

| Scenario | Biological background | Named override |
|---|---|---|
| `normal_human_plasma` | human plasma-like medium | none |
| `mouse_plasma` | mouse plasma/interstitial-fluid composition | none |
| `high_glucose` | same human plasma-like background | glucose increased |
| `low_glucose` | same human plasma-like background | glucose reduced |
| `high_lactate` | same human plasma-like background | lactate increased |
| `low_lactate` | same human plasma-like background | lactate reduced |
| `low_glutamine` | same human plasma-like background | glutamine reduced |

Challenge scenarios retain the same basal background and change only the named
component. They are modelling environments, not measured transporter fluxes.
Concentration-derived uptake caps remain explicit assumptions and are
intersected with the original GEM directionality.

Inspect composition and provenance fields:

```r
unique(medium_scenarios[, intersect(c(
  "medium_scenario_id",
  "medium_background_id",
  "composition_primary_reference_doi",
  "composition_validation_reference_doi",
  "challenge_reference_doi",
  "scenario_construction"
), colnames(medium_scenarios))])
```

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

Built-in and custom scenarios may be generated together. Full references and
interpretation rules are in [Medium scenarios and evidence](medium-presets.md).

## One-shot run

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
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
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

`upstream_workers` controls the GRN and other upstream parallel stages;
`layer2_workers` controls the LP-heavy metabolic stage. Layer 2 workers force
numerical libraries and HiGHS to one internal thread to avoid nested
oversubscription. On memory-limited systems, reduce `layer2_workers` before
changing the model definition.

Do not place `parallel` or `BPPARAM` inside `pando_infer_args`; the workflow owns
Stage 1 parallelism. The following retired condition-GRN controls are rejected:

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

For every retained broad cell type with at least two eligible conditions:

1. discover candidate TF–peak–target edges on all eligible cells of the cell
   type;
2. repeat candidate discovery separately in each condition;
3. union complete `(target, TF, region)` triples without Cartesian
   recombination;
4. freeze the target-specific dictionary;
5. fit every condition with the same Gaussian identity model
   `target ~ TF:peak`, using shared preprocessing and `scale = FALSE`;
6. retain the complete coefficient and uncertainty table;
7. pass only estimable coefficients with BH `padj < 0.05` to the condition
   penalty.

The multi-condition penalty route does not apply an additional absolute-effect
or model-R² threshold after Pando's significance flag. `min_abs_estimate` and
`min_model_rsq` remain relevant to the direct standard-Pando fallback where its
legacy edge extraction is used.

A zero-variance, aliased, non-finite, or insufficient-residual-df coefficient
remains `NA` and is excluded. A non-significant edge remains in the complete
table and must not be interpreted as a biological zero. GLM P values are
conditional on the frozen candidate dictionary.

## No condition or one condition

When `condition_col = NULL`, the requested column is absent, or only one
non-missing level is observed, Stage 1 runs the original Pando Gaussian
interaction GRN independently for each retained broad cell type. No condition
coefficient or synthetic condition fit is created.

## Metacells and penalty handoff

Stage 2 builds one independent multimodal WNN graph per broad cell type. All
conditions of that type share the graph, and condition splits parent membership
after clustering so final metacells remain condition-pure.

For multiple conditions, RegCompass computes each retained paired-cell
contribution as

```text
penalty_effect × TF_RNA × peak_ATAC
```

using the coefficient from that cell's own condition, sums contributions by
target, and only then averages over exact SuperCell membership. It does not
recompute TF×ATAC from metacell averages or refit coefficients after
aggregation.

`regulatory_alpha = 1` and `gpr_and_method = "min"` remain canonical. The same
medium-specific GEM, reaction order, bounds, target direction, and `vmax` are
reused across conditions and metacells.

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

Several Stage 4/5 fields retain historical names containing `_oof`, `common`, or
`condition_unique`. In the current condition model:

```text
gene_projection_condition_full_oof   = primary fixed-dictionary projection
gene_projection_common_oof           = compatibility alias of primary
gene_projection_condition_unique_oof = zero compatibility matrix
penalty_condition_full_oof            = primary penalty
penalty_common_oof                    = compatibility alias of primary
penalty_condition_unique_increment    = zero compatibility matrix
```

These names do not imply OOF estimation or a shared-slope decomposition.

Use [Tutorial 2](tutorial-02-stepwise-audit.md) for restartable stages,
[Tutorial 4](tutorial-04-targeted-reaction-remapping.md) for targeted reaction
remapping, and [Tutorial 5](tutorial-05-condition-differential-analysis.md) for
condition statistics.
