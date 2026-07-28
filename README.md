# RegCompassR

RegCompassR 1.9.1 implements an RNA+ATAC metabolic workflow for paired
single-cell multiome data.

## Workflow

```text
cell type cells across conditions
→ one shared Pando TF–peak–GEM-gene design
→ independently estimated, directly comparable condition coefficients
→ supported metabolic target genes
→ complete-GPR core reactions
→ one ordered subsystem/cross-reference expansion pass
→ integrated RNA+ATAC reaction support
→ medium-constrained model with global FASTCORE completion
→ directional COMPASS-like LP scoring
→ annotated rankings and condition contrasts
```

Pando is fitted once per cell type across all conditions. Conditions use the
same complete TF–peak–target edge dictionary, edge eligibility mask, and pooled
final-predictor scale. Within each target, they also share a lambda path and
selected lambda, while their elastic-net coefficients are estimated
independently. This is designed to retain the
interpretation of separate condition fits without allowing edge definitions or
coefficient units to drift.

The candidate target genes are all GEM GPR genes present in the RNA assay. A
gene enters the Stage 3 supported set when at least one active
condition-level TF–peak–gene coefficient passes the effect-size and
target-model-R² filters. Positive and negative coefficients both count as
regulatory evidence. A reaction is a core only when one complete GPR branch is
contained in that supported gene set.

For the condition-dependent penalty modifier, RegCompass uses Pando's explicit
reference contrast
`Δβ = β_condition - β_reference` and reconstructs the exact standardized
predictor from metacell TF RNA and peak ATAC. It does not divide edge effects by
their absolute sum, so fitted effect amplitude is preserved.

Stage 3 expansion is fixed and executed exactly once:

```text
complete-GPR cores
→ all reactions in core-reaction subsystems
→ direct KEGG/Reactome reaction equivalents
→ direct master-Rhea reaction equivalents
```

For Layer 1 reaction support, genes joined by a GPR AND relationship are aggregated with one of the three COMPASS functions: `min`, `median`, or `mean`. RegCompass defaults to `min`, representing the limiting required subunit. Isozyme OR branches remain additive in the canonical workflow.

## Installation

### Default validated profile: Seurat v4

The canonical and release-validation environment remains the pinned Seurat v4 stack. This is the default profile for reproducing RegCompass analyses and for creating new input objects.

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github(
  "1667857557/SuperCell_Seurat_V4@supercell-2.0"
)
remotes::install_github(
  "1667857557/Pando_regcompass@f98923101b00a138652479bddee44ab0b24b07b6"
)
remotes::install_github("1667857557/Regcompass")
```

### Optional compatible profile: Seurat v5

RegCompass also accepts SeuratObject/Seurat 5.x with Signac 1.x from version 1.12.0 onward. For the closest behavior to the default profile, create v3-style assays while running Seurat v5:

```r
options(Seurat.object.assay.version = "v3")
```

An existing v5 `Assay5` is supported when its `counts.*` and optional `data.*` layers can be joined without ambiguity. The canonical Stage 1 and Stage 2 functions join those layers in a working copy and record the operation in `object@misc$regcompass_seurat_compatibility`; the caller's original object is not rewritten. Signac 2.x `ChromatinAssay5` is not yet supported.

See [Seurat v4/v5 compatibility](docs/seurat-compatibility.md) for the version matrix, layer policy, provenance, and migration checks.

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

  # Stage 1
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      method = "shared_design_independent",
      candidate_screen = "condition_union",
      tf_cor = 0.1,
      peak_cor = 0.01,
      alpha = 0.5,
      condition_mix = 1,
      condition_weight = "equal",
      reference_condition = "Control",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE,
      parallel = FALSE
    )
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

Public runner arguments are ordered by processing sequence: shared model inputs → Stage 1 Pando → Stage 2 metacells → Stage 3 meta-modules → Stage 4 Layer 1 → Stage 5 Layer 2 → execution controls.

When `pfm` is omitted, RegCompass internally performs the equivalent of:

```r
data("motifs", package = "Pando")
pfm <- motifs
```

A user-supplied `pfm` still overrides the default.

The condition-comparable Pando defaults are
`method = "shared_design_independent"`,
`candidate_screen = "condition_union"`, `condition_mix = 1`,
`condition_weight = "equal"`, and `scale = TRUE`. RegCompass rejects
incompatible overrides because they would change the meaning or units of
condition contrasts. `min_abs_estimate = 0` and `min_model_rsq = 0.1` control
the downstream active-edge filter. Adjusted P values are not produced by this
regularized solver; the legacy `padj_threshold` and `require_padj` fields are
retained only for call compatibility.

The canonical metacell geometry defaults are RNA `pca` dimensions `1:30`, ATAC `lsi` dimensions `2:30`, and `seed = 12345L`. For ordered condition strata, the internal SuperCell2 seed is `seed + stratum_index - 1`. Changing cells, assay matrices, reductions, dimensions, seed, gamma, or metacell thresholds requires `overwrite = TRUE` to rebuild Stage 2 checkpoints.

`gpr_and_method` accepts only `"min"`, `"median"`, or `"mean"`; omitting it uses `"min"`.

`completion_time_limit` applies only while FASTCORE constructs the medium-specific union GEM. Directional scoring LPs run without a time-limit parameter.

Unless `pando_args$pando_initiate_args$regions` is supplied, RegCompass uses species-specific Pando region defaults:

```text
human: phastConsElements20Mammals.UCSC.hg38 ∪ SCREEN.ccRE.UCSC.hg38
mouse: phastConsElements20Mammals.UCSC.hg38 only
```

An explicit region object overrides either default.

## Medium presets

`rc_make_medium_scenarios()` supports physiological, culture-medium, nutrient-sensitivity, technical, and custom scenarios:

`physiologic`, `normal_human_plasma`, `mouse_plasma`, `rpmi1640`, `dmem_high_glucose`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, `minimal`, `compass_model_bounds`, `permissive_all_exchange`, and `custom`.

See [Predefined extracellular medium scenarios](docs/medium-presets.md) for species restrictions, assumptions, and custom-medium examples.

## Inspectable stages

- `rc_regcompass_step_grn()`: fit one shared-design, condition-comparable Pando model per cell type for GEM target genes.
- `rc_regcompass_step_metacells()`: construct condition-level multimodal metacells from explicit RNA/ATAC reductions, dimensions, and a reproducible seed.
- `rc_regcompass_step_meta_modules()`: summarize significant metabolic targets, map complete-GPR cores, and perform one fixed ordered annotation expansion pass.
- `rc_regcompass_step_layer1()`: calculate integrated RNA+ATAC reaction support with COMPASS-compatible GPR-AND aggregation.
- `rc_regcompass_step_layer2()`: build the medium-constrained model and run directional LP scoring.
- `rc_regcompass_step_results()`: assemble rankings, annotations, provenance, and contrasts.
- `rc_regcompass_step_target_union()`: remap selected core genes or reactions and score directly linked targets in the cached Stage 5 model.

## Main outputs

```r
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
result$merged_grn_meta_modules$merged_reaction_membership
result$microcompass$model_cache_summary
```

## Tutorials

- [Level 1: minimal one-shot run](docs/tutorial-01-quick-start.md)
- [Level 2: stepwise run](docs/tutorial-02-stepwise-audit.md)
- [Level 3: restart, sensitivity, and diagnostics](docs/tutorial-03-advanced-restart.md)
- [Level 4: targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Level 5: condition differential analysis](docs/tutorial-05-condition-differential-analysis.md)
- [Condition-comparable Pando contract and penalty projection](docs/condition_multitask_grn.md)
- [Seurat v4/v5 compatibility](docs/seurat-compatibility.md)
- [Metacell reduction and seed selection](docs/metacell-reduction-selection.md)
- [Medium presets](docs/medium-presets.md)
- [Workflow and mathematical interpretation](docs/workflow.md)
- [Stage input-output contracts](docs/stage-interface-contracts.md)
