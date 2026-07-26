# Tutorial Level 1: minimal one-shot run

Use this tutorial for paired single-cell RNA+ATAC data and RegCompassR 1.8.8.

## Workflow

```text
one shared Pando TF–peak–target background per cell type
→ global GRN backbone + condition deviations
→ stability-selected condition sub-GRNs
→ condition-specific metabolic target genes
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

The object must contain paired RNA and ATAC measurements, `condition_col`, `celltype_col`, and optionally a biological `sample_col`. RegCompass normalises RNA globally and ATAC with one TF-IDF reference per cell type across conditions.

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
  sample_col = "sample_id",

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
    n_stability = 50,
    stability_fraction = 0.8,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    seed = 12345L
  ),

  # Stage 2: condition-level metacells
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

## Stage 1 parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `min_cells` | `20L` | Minimum cells in every condition of a cell type. |
| `alpha` | `0.5` | Elastic-net mixing value. It must be below one so the symmetric deviation solution has a ridge component. |
| `global_penalty_factor` | `1` | Penalty factor for the shared backbone. |
| `deviation_penalty_factor` | `2` | Stronger shrinkage for condition deviations. |
| `n_stability` | `25L` | Repeated stratified subsamples. Use a larger value for final analyses. |
| `min_selection_frequency` | `0.7` | Minimum fraction of successful subsamples selecting the condition edge. |
| `min_sign_stability` | `0.8` | Minimum conditional sign agreement. |
| `candidate_screen_threshold` | `0` | Default retains the complete structural Pando candidate universe. |
| `max_edges_per_target` | `Inf` | Default does not truncate the shared candidate universe. |

The multitask model does not assign classical adjusted p-values to regularised coefficients. `padj` is `NA`; edge support is defined by stability selection and target-model cross-validated reliability.

## Core reaction rule

For condition target set `G_c`, reaction `r` becomes a core only when one complete GPR branch is present:

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
