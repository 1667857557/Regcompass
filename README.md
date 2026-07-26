# RegCompassR

RegCompassR 1.8.8 implements a condition-comparable RNA+ATAC regulatory–metabolic workflow for paired single-cell multiome data.

## Workflow

```text
all conditions within one cell type
→ one shared Pando structural TF–peak–GEM-gene candidate universe
→ condition-balanced global-plus-deviation elastic-net model
→ stability-selected condition sub-GRNs
→ condition-specific regulated metabolic target genes
→ condition-specific complete-GPR core reactions
→ one ordered subsystem/cross-reference expansion pass
→ merged reaction catalogue across all conditions and cell types
→ integrated RNA+ATAC reaction support
→ one shared medium-specific union GEM with global FASTCORE completion
→ directional COMPASS-like LP scoring on the identical structure
→ annotated rankings and condition contrasts
```

The default Stage 1 mode is `grn_mode = "multitask_shared_backbone"`. Pando is called once per cell type to construct a common structural edge universe

```text
U_m = {(TF, regulatory peak, metabolic target gene)}
```

before coefficients are estimated. Every condition therefore uses identical candidate columns. RegCompass then estimates

```text
theta[e, c] = beta[e] + delta[e, c]
```

where `beta` is the cell-type global backbone and `delta` is the condition deviation. The reported decomposition is reference-free:

```text
beta[e]     = mean over conditions of theta[e, c]
delta[e, c] = theta[e, c] - beta[e]
```

Positive and negative stable edges both identify a condition-regulated metabolic gene. A reaction is a condition core only when at least one complete GPR AND branch is contained in that condition/cell-type target set.

The former independent `condition × cell type` Pando workflow remains available as:

```r
grn_mode = "legacy_condition_pando"
```

## Biological meta-module expansion

Stage 3 expansion is fixed and executed exactly once:

```text
complete-GPR cores
→ all reactions in core-reaction subsystems
→ direct KEGG/Reactome reaction equivalents
→ direct master-Rhea reaction equivalents
```

Stage 3 does not run FASTCORE and does not construct a GEM. For each medium, Stage 5 merges all condition/cell-type reaction memberships and constructs one shared union GEM. The same stoichiometric matrix, bounds and target catalogue are reused for every condition and metacell.

For Layer 1 reaction support, genes joined by a GPR AND relationship are aggregated with `min`, `median`, or `mean`. RegCompass defaults to `min`, representing the limiting required subunit. Isozyme OR branches remain additive.

## Installation

### Default validated profile: Seurat v4

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github(
  "1667857557/SuperCell_Seurat_V4@supercell-2.0"
)
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

RegCompass 1.8.8 requires the Pando fork providing `prepare_grn_design()`.

### Optional compatible profile: Seurat v5

RegCompass also accepts SeuratObject/Seurat 5.x with Signac 1.x from version 1.12.0 onward. For behavior closest to the default profile:

```r
options(Seurat.object.assay.version = "v3")
```

An existing v5 `Assay5` is supported when its `counts.*` and optional `data.*` layers can be joined without ambiguity. The caller's original object is not rewritten. Signac 2.x `ChromatinAssay5` is not supported.

See [Seurat v4/v5 compatibility](docs/seurat-compatibility.md).

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
  sample_col = "sample_id",

  # Stage 1: shared candidate universe and joint condition model
  grn_mode = "multitask_shared_backbone",
  pando_args = list(
    min_cells = 300,
    pando_design_args = list(
      screen_method = "structural",
      min_tf_detection = 0.01,
      min_peak_detection = 0.01,
      min_target_detection = 0.01
    )
  ),
  multitask_args = list(
    alpha = 0.5,
    global_penalty_factor = 1,
    deviation_penalty_factor = 2,
    lambda_rule = "lambda.1se",
    nfolds = 5,
    n_stability = 30,
    stability_fraction = 0.8,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    min_abs_effect = 0,
    min_cv_rsq = 0.1,
    seed = 12345L
  ),

  # Stage 2
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

  # Stage 4
  layer1_args = list(
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),

  # Stage 5
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

When `pfm` is omitted, RegCompass loads:

```r
data("motifs", package = "Pando")
pfm <- motifs
```

Unless `pando_args$pando_initiate_args$regions` is supplied, the species-specific region defaults are:

```text
human: phastConsElements20Mammals.UCSC.hg38 ∪ SCREEN.ccRE.UCSC.hg38
mouse: phastConsElements20Mammals.UCSC.hg38 only
```

### Shared-backbone model controls

- `alpha` must be below 1 so the ridge component makes the redundant global/deviation design uniquely estimable.
- `deviation_penalty_factor` must exceed `global_penalty_factor`, shrinking weak condition-specific effects toward the common backbone.
- observation weights give every condition equal total loss weight.
- folds are stratified by condition.
- stability resampling is stratified by condition or by condition–sample intercept stratum.
- `selection_frequency` and `sign_stability` are stability diagnostics, not classical p-values.

The default structural Pando design does not use pooled target-RNA significance. This prevents opposite condition effects from cancelling before the joint model is fitted. Optional `screen_method = "union_within_group_correlation"` retains a shared universe but requires a candidate to pass within at least one condition.

### Layer 1 regulatory projection

For multitask edges, RegCompass projects ATAC-only deviations using the effective condition coefficient, stability weight, a cell-type-shared TF reference and the interaction scale used during fitting. TF edges sharing one peak are signed-summed before the peak is projected, preventing repeated use of the same accessibility measurement. One denominator is shared across conditions for each target and cell type, preserving relative regulatory strength.

The bounded RNA integration remains:

```text
C_multiome = C_RNA * 2^(alpha * R_ATAC) /
             (1 - C_RNA + C_RNA * 2^(alpha * R_ATAC))
```

Therefore `R_ATAC = 0` returns exactly the RNA-only support.

## Legacy Stage 1 example

```r
result_legacy <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_legacy",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "Group",
  celltype_col = "cell_type",
  grn_mode = "legacy_condition_pando",
  pando_args = list(
    min_cells = 300,
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = TRUE,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    )
  )
)
```

## Main Stage 1 outputs

```r
result$grn$tf_peak_gene_candidates
result$grn$tf_peak_gene_global
result$grn$tf_peak_gene_condition_all
result$grn$tf_peak_gene_significant
result$grn$condition_target_genes
result$grn$celltype_fit_status
```

The condition edge table includes:

```text
global_estimate
condition_deviation
effective_estimate
selection_frequency
sign_stability
stability_weight
active_edge
effective_direction
sign_flip_flag
design_id
```

## Main metabolic outputs

```r
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
result$merged_grn_meta_modules$merged_reaction_membership
result$microcompass$model_cache_summary
```

## Medium presets

`rc_make_medium_scenarios()` supports `physiologic`, `normal_human_plasma`, `mouse_plasma`, `rpmi1640`, `dmem_high_glucose`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, `minimal`, `compass_model_bounds`, `permissive_all_exchange`, and `custom`.

`completion_time_limit` applies only to FASTCORE construction of the medium-specific union GEM. Directional scoring LPs have no scoring time-limit argument.

## Inspectable stages

- `rc_regcompass_step_grn()`: construct shared cell-type Pando designs and fit condition-comparable multitask sub-GRNs.
- `rc_regcompass_step_metacells()`: construct condition-level multimodal metacells.
- `rc_regcompass_step_meta_modules()`: map active condition target genes to complete-GPR cores and annotation-expanded modules.
- `rc_regcompass_step_layer1()`: calculate integrated RNA+ATAC reaction support.
- `rc_regcompass_step_layer2()`: construct the shared medium-specific union GEM and run directional LP scoring.
- `rc_regcompass_step_results()`: assemble rankings, annotations, provenance and contrasts.
- `rc_regcompass_step_target_union()`: remap selected genes or reactions and score linked targets in the cached Stage 5 model.

## Tutorials

- [Shared-backbone GRN mathematics and object contracts](docs/multitask-shared-grn.md)
- [Level 1: minimal one-shot run](docs/tutorial-01-quick-start.md)
- [Level 2: stepwise run](docs/tutorial-02-stepwise-audit.md)
- [Level 3: restart, sensitivity, and diagnostics](docs/tutorial-03-advanced-restart.md)
- [Level 4: targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Level 5: condition differential analysis](docs/tutorial-05-condition-differential-analysis.md)
- [Seurat v4/v5 compatibility](docs/seurat-compatibility.md)
- [Metacell reduction and seed selection](docs/metacell-reduction-selection.md)
- [Medium presets](docs/medium-presets.md)
- [Workflow and mathematical interpretation](docs/workflow.md)
- [Stage input-output contracts](docs/stage-interface-contracts.md)
