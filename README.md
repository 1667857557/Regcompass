# RegCompassR

RegCompassR 1.8.4 implements a GRN-first RNA+ATAC metabolic workflow for paired single-cell multiome data.

## Canonical architecture

```text
condition × cell type single cells
→ one Pando GRN per group
→ condition-level multimodal metacells
→ metabolic-gene GRN components
→ complete-GPR core reactions
→ core-subsystem + KEGG/Reactome + master-Rhea expansion
→ biological meta-modules
→ deduplicated merged meta-module catalogue
→ integrated RNA+ATAC reaction support
→ one medium-specific union GEM
→ one global FASTCORE completion
→ directional two-step COMPASS-like LP scoring
```

### Terminology

- **Meta-module**: a biological reaction set defined from the GRN, complete GPRs, subsystem membership, and direct reaction cross-references.
- **Merged meta-module catalogue**: the reaction-ID deduplication of all biological meta-modules. It is not flux-completed and is not a GEM.
- **Union GEM**: the medium-constrained Stage 5 model created from the merged catalogue plus global FASTCORE support. Only this model is called a union GEM.

FASTCORE is applied once for each medium-specific union GEM. Biological meta-modules and their merged catalogue are never FASTCORE-completed independently.

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

## Inspectable stages

- `rc_regcompass_step_grn()`: condition-by-cell-type Pando GRNs.
- `rc_regcompass_step_metacells()`: condition-level, cell-type-guided SuperCell2 metacells.
- `rc_regcompass_step_meta_modules()`: complete-GPR cores and biological reaction expansion; no FASTCORE and no GEM construction.
- `rc_regcompass_step_layer1()`: integrated RNA+ATAC reaction support.
- `rc_regcompass_step_layer2()`: medium-specific union-GEM construction, global FASTCORE, and directional LP scoring.
- `rc_regcompass_step_results()`: rankings, evidence provenance, annotations, and condition contrasts.
- `rc_regcompass_step_target_union()`: second-pass scoring of directly KEGG/Reactome/master-Rhea-linked non-core reactions in the exact cached final medium-specific union GEMs.

## Key object fields

```r
step3$condition_modules
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership

step5$model_cache_summary
result$merged_grn_meta_modules
result$reaction_ranking
result$condition_contrast
```

## Structural interpretation

For one medium scenario, all conditions and metacells use the same union GEM, reaction bounds, target-flux fraction, and direction-specific `vmax`. Condition differences arise from the multiome penalty matrix.

Different medium scenarios may produce different global FASTCORE support sets and therefore different union GEM structures. Cross-medium comparisons must be interpreted as different structural contexts.

The optional second-pass remapping workflow uses selected original core reactions as database-mapping anchors. It reuses the exact Stage 5 union-GEM files, medium bounds, and reaction structure; it does not rebuild a model or rerun FASTCORE.

## Tutorials

- [Level 1: minimal one-shot run](docs/tutorial-01-quick-start.md)
- [Level 2: true stepwise run with audit gates](docs/tutorial-02-stepwise-audit.md)
- [Level 3: restart, sensitivity, and diagnostics](docs/tutorial-03-advanced-restart.md)
- [Level 4: targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Level 5: condition differential analysis](docs/tutorial-05-condition-differential-analysis.md)
- [Workflow and mathematical interpretation](docs/workflow.md)
- [Stage input-output contracts](docs/stage-interface-contracts.md)
