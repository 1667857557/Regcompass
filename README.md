# RegCompassR

RegCompassR 1.8.4 implements a GRN-first RNA+ATAC metabolic workflow for paired single-cell multiome data.

## Workflow

```text
condition × cell type cells
→ group-specific Pando GRNs and condition-level metacells
→ complete-GPR reaction meta-modules
→ integrated RNA+ATAC reaction support
→ one model per medium with one global FASTCORE completion
→ directional COMPASS-like LP scoring
→ annotated rankings and condition contrasts
```

## Installation

```r
remotes::install_github("1667857557/Regcompass")
```

Pando must be installed from the RegCompass-compatible fork:

```r
remotes::install_github("1667857557/Pando_regcompass")
```

## Minimal complete run

```r
library(RegCompassR)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motif2tf, package = "Pando")

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
  pfm = motif2tf,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  gem = gem,
  species = "human",
  medium_scenarios = medium_scenarios,
  sample_col = NULL,
  condition_col = "Group",
  celltype_col = "cell_type",
  model_mode = "meta_module_gem",
  metacell_args = list(
    gamma = 30,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  pando_args = list(
    min_cells = 100,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    )
  ),
  layer1_args = list(
    top_k_neighbors = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    max_targets_per_tf = 200,
    expansion_mode = "ordered_once",
    regulatory_alpha = 1,
    tau = 0.20
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 600,
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

## Medium presets

`rc_make_medium_scenarios()` supports physiological, culture-medium, nutrient-sensitivity, technical, and custom scenarios:

`physiologic`, `normal_human_plasma`, `mouse_plasma`, `rpmi1640`, `dmem_high_glucose`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, `minimal`, `compass_model_bounds`, `permissive_all_exchange`, and `custom`.

See [Predefined extracellular medium scenarios](docs/medium-presets.md) for species restrictions, assumptions, and custom-medium examples.

## Inspectable stages

- `rc_regcompass_step_grn()`: infer condition-by-cell-type Pando GRNs.
- `rc_regcompass_step_metacells()`: construct condition-level multimodal metacells.
- `rc_regcompass_step_meta_modules()`: map GRN components to complete-GPR reaction modules.
- `rc_regcompass_step_layer1()`: calculate integrated RNA+ATAC reaction support.
- `rc_regcompass_step_layer2()`: build medium-specific models and run directional LP scoring.
- `rc_regcompass_step_results()`: assemble rankings, annotations, provenance, and contrasts.
- `rc_regcompass_step_target_union()`: remap selected core genes or reactions and score directly linked targets in the cached Stage 5 models.

## Main outputs

```r
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
result$merged_grn_meta_modules$merged_reaction_membership
result$microcompass$model_cache_summary
```

Within one medium scenario, all conditions use the same structural model and reaction bounds; condition differences arise from the multiome penalty matrix. Different media may produce different FASTCORE support sets and should be interpreted as different structural contexts.

## Tutorials

- [Level 1: minimal one-shot run](docs/tutorial-01-quick-start.md)
- [Level 2: stepwise run](docs/tutorial-02-stepwise-audit.md)
- [Level 3: restart, sensitivity, and diagnostics](docs/tutorial-03-advanced-restart.md)
- [Level 4: targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Level 5: condition differential analysis](docs/tutorial-05-condition-differential-analysis.md)
- [Medium presets](docs/medium-presets.md)
- [Workflow and mathematical interpretation](docs/workflow.md)
- [Stage input-output contracts](docs/stage-interface-contracts.md)
