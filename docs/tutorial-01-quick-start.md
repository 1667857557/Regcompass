# Tutorial Level 1: minimal one-shot run

Use this tutorial for paired single-cell RNA+ATAC data and RegCompassR 1.8.8.

## Workflow

```text
one validated Pando TF–peak–target background per cell type
→ global GRN backbone + symmetric condition deviations
→ full-size condition-stratified bootstrap stability
→ condition-specific sub-GRNs and metabolic target genes
→ complete-GPR condition core reactions
→ one ordered subsystem/cross-reference expansion pass
→ one shared medium-specific union GEM
→ RNA+ATAC penalties and directional LP scoring
```

## Prepare the model

```r
library(RegCompassR)
library(Seurat)
library(Signac)
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
```

The object must contain paired RNA and ATAC measurements and complete metadata columns supplied as `condition_col` and `celltype_col`. The canonical workflow does not accept a biological-sample column. RNA is normalised globally, and ATAC uses one TF-IDF reference per cell type across conditions.

When `pfm` is omitted, RegCompass loads `data("motifs", package = "Pando")`. Default regulatory regions are:

```text
human: phastConsElements20Mammals.UCSC.hg38 ∪ SCREEN.ccRE.UCSC.hg38
mouse: phastConsElements20Mammals.UCSC.hg38
```

## Run the complete workflow

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,

  condition_col = "Group",
  celltype_col = "cell_type",

  # Stage 1: shared candidate background and condition sub-GRNs
  grn_mode = "multitask_shared_backbone",
  pando_args = list(
    min_cells = 100,
    pando_design_args = list(
      peak_to_gene_method = "Signac",
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
    n_bootstrap = 100,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    seed = 12345L
  ),

  # Stage 2: condition-only, cell-type-label-guided metacells
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30,
    seed = 12345L,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  ),

  # Stage 4: integrated evidence and GPR aggregation
  layer1_args = list(
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),

  # Stage 5: one shared model per medium
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
  layer2_workers = 30,
  progress = TRUE
)
```

The package default is `n_bootstrap = 50L`; `100` is shown for a final analysis with a more stable empirical selection-frequency estimate.

## Stage 1 parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `min_cells` | `20L` | Minimum cells in every condition of a cell type. |
| `alpha` | `0.5` | Elastic-net mixing value; it must remain below one so the symmetric deviation solution contains a ridge component. |
| `global_penalty_factor` | `1` | Penalty factor for the shared backbone. |
| `deviation_penalty_factor` | `2` | Stronger default shrinkage for condition deviations. |
| `nfolds` | `5L` | Maximum number of condition-stratified cell-level CV folds. |
| `n_bootstrap` | `50L` | Full-size nonparametric bootstrap replicates, sampled with replacement within each condition. |
| `min_selection_frequency` | `0.7` | Minimum fraction of successful bootstrap fits selecting the condition edge. |
| `min_sign_stability` | `0.8` | Minimum conditional sign agreement among selected bootstrap fits. |
| `candidate_screen_threshold` | `0` | Retains the complete structural Pando candidate universe by default. |
| `max_edges_per_target` | `Inf` | Does not truncate the shared candidate universe by default. |

For every bootstrap replicate, each condition contributes exactly its original cell count, sampled with replacement. The resampled target and TF-by-ATAC predictors are re-centred within condition before refitting at the full-data selected lambda.

The multitask model does not assign classical adjusted p-values to regularised coefficients. `padj` is `NA`; edge support is defined by bootstrap selection frequency, conditional sign stability, the full-data absolute effect, and target-model cross-validated reliability.

## Core reaction rule

For condition target set `G_c`, reaction `r` becomes core only when one complete GPR branch is present:

\[
Core_{r,c}=1\iff\exists k:B_{r,k}\subseteq G_c.
\]

Positive and negative active regulatory edges both allow a target gene to enter `G_c`.

## Inspect outputs

```r
result$grn$tf_peak_gene_candidates
result$grn$tf_peak_gene_global
result$grn$tf_peak_gene_condition_all
result$grn$tf_peak_gene_significant
result$grn$condition_target_genes
result$grn$stability_diagnostics
result$grn$group_status

result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$merged_grn_meta_modules$merged_core_reactions
result$merged_grn_meta_modules$merged_reaction_membership

result$reaction_ranking
result$condition_contrast
result$microcompass$model_cache_summary
```

The merged catalogue is a reaction union, not a GEM. Stage 5 constructs one medium-specific union GEM and reuses its exact stoichiometry and bounds for every condition and metacell.

See [multitask GRN mathematics and stage contracts](multitask-shared-grn.md).
