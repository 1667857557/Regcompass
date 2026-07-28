# Public functions and API contract in RegCompassR 1.9.3

RegCompass exposes a one-shot workflow and six restartable stages. Pando 1.2.1
is the sole condition-GRN estimator; RegCompass consumes its versioned
`ConditionGRNFit` and does not refit condition coefficient matrices.

## Complete workflows

| Function | Purpose |
|---|---|
| `rc_run_regcompass_one_shot()` | Load species-aware GEM/medium defaults and run all stages. |
| `rc_run_regcompass()` | Run all stages with an explicit GEM, medium, and worker budgets. |

The complete workflow accepts `pando_args`, `metacell_args`,
`meta_module_args`, `layer1_args`, `layer2_args`, `upstream_workers`, and
`layer2_workers`.

## Restartable stages

| Stage | Function | Main output |
|---:|---|---|
| 1 | `rc_regcompass_step_grn()` | Pando `ConditionGRNFit` objects plus condition coefficient/effect tables |
| 2 | `rc_regcompass_step_metacells()` | condition-level multimodal metacells and cache provenance |
| 3 | `rc_regcompass_step_meta_modules()` | complete-GPR cores and merged biological reaction catalogue |
| 4 | `rc_regcompass_step_layer1()` | RNA-only and RNA+ATAC reaction support |
| 5 | `rc_regcompass_step_layer2()` | shared medium-specific model and directional LP scores |
| 6 | `rc_regcompass_step_results()` | annotated rankings, contrasts, evidence, and provenance |

## Stage 1 Pando contract

```r
pando_args = list(
  min_cells = 20L,
  min_abs_estimate = 0,
  min_model_rsq = 0.1,
  pando_infer_args = list(
    method = "shared_design_independent",
    candidate_screen = "motif_domain",
    condition_mix = 1,
    condition_weight = "equal",
    reference_condition = "Control",
    scale = TRUE
  )
)
```

Within each cell type, Pando shares:

- one TF–peak–target edge dictionary;
- pooled final `TF RNA × peak ATAC` predictor transforms;
- one pooled target transform;
- one target-specific lambda path and selected lambda.

Condition coefficient columns are estimated independently at the selected
lambda. `candidate_screen = "motif_domain"` is the interaction-safe default.
The optional `condition_union` and `pooled` modes impose marginal-correlation
screens and are intended for explicit sensitivity analyses.

### Comparison support

Pando exports both `eligibility_mask` and `comparison_mask`:

```text
comparison_mask[e, c] =
  eligibility_mask[e, c] && eligibility_mask[e, reference]
```

RegCompass requires the explicit Pando 1.2.1 mask. All effect rows expose
`comparable_to_reference`; only comparable rows can enter
`tf_peak_gene_condition_effect` and the Layer 1 regulatory projection.

### Pando argument routing

| Nested bundle | Target function |
|---|---|
| `pando_initiate_args` | `Pando::initiate_grn()` |
| `pando_motif_args` | `Pando::find_motifs()` |
| `pando_infer_args` | `Pando::infer_condition_grn()` |

RegCompass rejects nested overrides of:

```text
object, peak_assay, rna_assay, pfm, genome,
cell_type_col, condition_col, genes, network_name,
min_cells_per_condition, on_small_condition, BPPARAM
```

It also rejects `aggregate_rna_col` and `aggregate_peaks_col` because canonical
Stage 1 fits paired single cells and Stage 2 owns aggregation.

### Stage 1 parallel execution

For `rc_regcompass_step_grn()`:

| `parallel` | `BPPARAM` | Pando route |
|---|---|---|
| `FALSE` | any | serial; a supplied backend is ignored with a warning |
| `TRUE` | `BiocParallelParam` | Pando receives that backend |
| `TRUE` | `NULL`/`FALSE` | Pando native map backend |

`BPPARAM = TRUE` is invalid. The resolved route is stored in
`step1$params$pando_parallel`. Complete workflows use `upstream_workers` and a
stage-scoped backend; users should not set `pando_infer_args$parallel`.

### Stage 1 outputs

```r
step1$grn_result$condition_grn_fits
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition_all
step1$grn_result$tf_peak_gene_condition
step1$grn_result$tf_peak_gene_condition_effect_all
step1$grn_result$tf_peak_gene_condition_effect
step1$grn_result$normalization_policy
step1$params$pando_parallel
```

Compatibility aliases `sample_status`, `tf_peak_gene_all`, and
`tf_peak_gene_significant` remain available but should not be used by new code.

## Regulatory regions

Human Stage 1 defaults to the union of Pando's hg38 phastCons and SCREEN ccRE
objects. Pando does not bundle a mouse-coordinate equivalent. Mouse runs must
supply a build-matched `GRanges` through
`pando_args$pando_initiate_args$regions`; RegCompass stops rather than applying
hg38 regions to mouse peaks.

## Stage 2 metacells

```r
metacell_args = list(
  rna_reduction = "pca",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30L,
  seed = 12345L,
  min_cells_per_stratum = 100L,
  min_metacell_size = 20L,
  min_metacells_per_stratum = 2L
)
```

Reduction names, dimensions, embedding fingerprints, ordered cells, assays,
seed, gamma, and thresholds are part of the cache contract.

## Stages 3–5

Stage 3 maps active condition coefficients to supported metabolic genes,
complete-GPR cores, and one ordered expansion pass:

```text
core subsystems → direct KEGG/Reactome equivalents → direct master-Rhea equivalents
```

Stage 4 reconstructs Pando's stored standardized interaction predictors and
projects only comparable condition-versus-reference coefficients. GPR AND
accepts `"min"`, `"median"`, or `"mean"`; the default is `"min"`. Isozyme OR
branches remain additive.

With `model_mode = "meta_module_gem"`, Stage 5 applies each medium, performs one
global FASTCORE completion, caches one shared model, and reuses it for every
condition and metacell in that medium.

## Other public analysis functions

| Function | Purpose |
|---|---|
| `rc_prepare_gem()` | Load/validate a pinned Human-GEM or Mouse-GEM. |
| `rc_make_medium_scenarios()` | Construct condition-invariant medium constraints. |
| `rc_build_reaction_annotations()` | Build reaction names, formulas, GPRs, and database cross-references. |
| `rc_attach_reaction_annotations()` | Attach the current annotation contract to an existing result. |
| `rc_select_gene_reactions()` | Select scored reactions by metabolic gene. |
| `rc_test_condition_reactions()` | Compare fixed reaction-direction targets across conditions. |
| `rc_report_condition_directions()` | Retain forward/reverse targets and derive support summaries. |
| `rc_plot_condition_reaction()` | Plot one reaction-direction target across conditions. |
| `rc_plot_condition_gene_reactions()` | Select and plot reaction targets associated with genes. |
| `rc_regcompass_step_target_union()` | Score directly database-linked non-core targets in cached Stage 5 models. |

See [Tutorial 1](tutorial-01-quick-start.md),
[Tutorial 2](tutorial-02-stepwise-audit.md), and the
[Pando contract](condition-comparable-grn.md).
