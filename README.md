# RegCompassR

RegCompassR 1.8.4 implements an RNA+ATAC metabolic workflow for paired single-cell multiome data.

## Workflow

```text
condition × cell type cells
→ Pando TF–peak–GEM-gene models
→ significantly supported metabolic target genes
→ complete-GPR core reactions
→ one ordered subsystem/cross-reference expansion pass
→ integrated RNA+ATAC reaction support
→ medium-constrained model with global FASTCORE completion
→ directional COMPASS-like LP scoring
→ annotated rankings and condition contrasts
```

Pando is fitted separately for each `condition × cell type`. The candidate target genes are all GEM GPR genes present in the RNA assay. A gene enters the Stage 3 supported set when at least one TF–peak–gene coefficient passes the configured adjusted-P-value, effect-size, and target-model-R² filters. Positive and negative coefficients both count as regulatory evidence. A reaction is a core only when one complete GPR branch is contained in that supported gene set.

Stage 3 expansion is fixed and executed exactly once:

```text
complete-GPR cores
→ all reactions in core-reaction subsystems
→ direct KEGG/Reactome reaction equivalents
→ direct master-Rhea reaction equivalents
```

There is no `expansion_mode`, fixed-point recursion, `max_iterations`, one-hop reaction expansion, or stoichiometric-neighbour expansion.

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
remotes::install_github("1667857557/Pando_regcompass")
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

The canonical Pando evidence defaults are `padj_threshold = 0.05`, `min_abs_estimate = 0`, `min_model_rsq = 0.1`, and `require_padj = TRUE`. All conditions must pass for a TF–peak–target row to enter the supported metabolic-gene set.

The canonical metacell geometry defaults are RNA `pca` dimensions `1:30`, ATAC `lsi` dimensions `2:30`, and `seed = 12345L`. For ordered condition strata, the internal SuperCell2 seed is `seed + stratum_index - 1`. Changing cells, assay matrices, reductions, dimensions, seed, gamma, or metacell thresholds requires `overwrite = TRUE` to rebuild Stage 2 checkpoints.

`gpr_and_method` accepts only `"min"`, `"median"`, or `"mean"`; omitting it uses `"min"`. The retired Boltzmann soft-min and `tau` API are not supported.

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

- `rc_regcompass_step_grn()`: fit condition-by-cell-type Pando models for GEM target genes using Pando's bundled `motifs` and species-specific default regions.
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
- [Seurat v4/v5 compatibility](docs/seurat-compatibility.md)
- [Metacell reduction and seed selection](docs/metacell-reduction-selection.md)
- [Medium presets](docs/medium-presets.md)
- [Workflow and mathematical interpretation](docs/workflow.md)
- [Stage input-output contracts](docs/stage-interface-contracts.md)
