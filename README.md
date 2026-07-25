# RegCompassR

RegCompassR 1.8.4 implements an RNA+ATAC metabolic workflow for paired single-cell multiome data.

## Workflow

```text
condition × cell type cells
→ Pando TF–peak–Human-GEM-gene models
→ significantly supported metabolic target genes
→ complete-GPR core reactions
→ subsystem and reaction-equivalence biological expansion
→ integrated RNA+ATAC reaction support
→ medium-constrained model with global FASTCORE completion
→ directional COMPASS-like LP scoring
→ annotated rankings and condition contrasts
```

Pando is fitted separately for each `condition × cell type`. The candidate target genes are all Human-GEM GPR genes present in the RNA assay. A gene enters the Stage 3 supported set when at least one TF–peak–gene coefficient passes the configured adjusted-P-value, effect-size, and target-model-R² filters. Positive and negative coefficients both count as regulatory evidence. A reaction is a core only when one complete GPR branch is contained in that supported gene set.

## Installation

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

## Minimal complete run

```r
library(RegCompassR)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

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
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  gem = gem,
  species = "human",
  medium_scenarios = medium_scenarios,
  condition_col = "Group",
  celltype_col = "cell_type",
  model_mode = "meta_module_gem",
  pando_args = list(
    min_cells = 300,
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    )
  ),
  metacell_args = list(
    gamma = 30,
    min_cells_per_stratum = 300,
    min_metacell_size = 10
  ),
  meta_module_args = list(
    expansion_mode = "ordered_once"
  ),
  layer1_args = list(
    regulatory_alpha = 1,
    tau = 0.20
  ),
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

`completion_time_limit` applies only while FASTCORE constructs the medium-specific union GEM. Directional scoring LPs run without a time-limit parameter.

Unless `pando_args$pando_initiate_args$regions` is supplied, human analyses load the two Pando data objects below and pass their union to `Pando::initiate_grn()`:

```r
data("phastConsElements20Mammals.UCSC.hg38", package = "Pando")
data("SCREEN.ccRE.UCSC.hg38", package = "Pando")
regions <- union(
  phastConsElements20Mammals.UCSC.hg38,
  SCREEN.ccRE.UCSC.hg38
)
```

The bundled default is hg38-specific. Non-human analyses must supply an appropriate `pando_initiate_args$regions` object.

## Medium presets

`rc_make_medium_scenarios()` supports physiological, culture-medium, nutrient-sensitivity, technical, and custom scenarios:

`physiologic`, `normal_human_plasma`, `mouse_plasma`, `rpmi1640`, `dmem_high_glucose`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, `minimal`, `compass_model_bounds`, `permissive_all_exchange`, and `custom`.

See [Predefined extracellular medium scenarios](docs/medium-presets.md) for species restrictions, assumptions, and custom-medium examples.

## Inspectable stages

- `rc_regcompass_step_grn()`: fit condition-by-cell-type Pando models for Human-GEM target genes.
- `rc_regcompass_step_metacells()`: construct condition-level multimodal metacells.
- `rc_regcompass_step_meta_modules()`: summarize significant metabolic targets, map complete-GPR cores, and perform annotation expansion.
- `rc_regcompass_step_layer1()`: calculate integrated RNA+ATAC reaction support.
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
- [Medium presets](docs/medium-presets.md)
- [Workflow and mathematical interpretation](docs/workflow.md)
- [Stage input-output contracts](docs/stage-interface-contracts.md)
