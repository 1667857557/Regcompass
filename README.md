# RegCompassR

RegCompassR 1.8.10 implements a shared-background regulatory–metabolic workflow for paired single-cell RNA+ATAC data.

## Canonical architecture

```text
all conditions within one cell type
→ one validated GREAT-domain Pando structural TF–peak–metabolic-gene universe
→ condition-aware TF×peak/target observability filter
→ one shared model edge universe across conditions
→ condition-balanced direct condition-theta elastic net
→ condition-stratified sample/donor cluster bootstrap when sample_col is valid
→ explicit warning and condition-stratified cell-bootstrap fallback otherwise
→ condition-specific active sub-GRNs and metabolic targets
→ complete-GPR condition core reactions
→ ordered subsystem / KEGG–Reactome / master-Rhea expansion
→ merged biological reaction catalogue
→ one shared medium-specific union GEM reused by every condition and metacell
→ RNA+ATAC penalties
→ directional COMPASS-like LP scoring
→ compact final analysis tables
```

For edge `e = (TF, peak, target)`, RegCompass estimates a condition-specific coefficient `theta[e,c]`. All conditions of one cell type use the same structural candidates, filtered model dictionary, edge scale, penalty structure, and candidate ordering.

## Installation

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github(
  "1667857557/SuperCell_Seurat_V4@c8b94949cd8a5ff7403f9f186c516f8efbac9b6f"
)
remotes::install_github(
  "1667857557/Pando_regcompass@6f42c8143bec6610b001e714a51627337f6d9ba9"
)
remotes::install_github("1667857557/Regcompass")
```

## Metadata contract

Required metadata:

```text
condition_col
celltype_col
```

Recommended metadata:

```text
sample_col
```

`sample_col` is used only for Stage 1 bootstrap stability. When valid, samples/donors are sampled with replacement separately within each condition and every selected sample contributes all of its cells. This is a cluster bootstrap, so bootstrap cell counts may vary with sample size.

When `sample_col = NULL`, or the named metadata column does not exist, RegCompass prints the exact fallback reason and uses condition-stratified cell resampling. An existing sample column with missing or empty IDs is an error. Conditions with fewer than two unique samples retain cluster bootstrap but emit a low-replication warning.

Stage 2 remains condition-only:

```text
strata_cols = condition_col
label = celltype_col
```

There is no artificial condition-pool metadata field. The sample column does not alter metacell strata or reconstruct separate metabolic models.

## Minimal run

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
  grn_mode = "multitask_shared_backbone",
  pando_args = list(
    min_cells = 100,
    pando_design_args = list(
      peak_to_gene_method = "GREAT",
      upstream = 100000,
      downstream = 0,
      extend = 1000000,
      only_tss = FALSE,
      min_tf_detection = 0,
      min_peak_detection = 0,
      min_target_detection = 0,
      max_edges_per_target = Inf
    )
  ),
  multitask_args = list(
    alpha = 0.5,
    global_penalty_factor = 1,
    deviation_penalty_factor = 1,
    lambda_rule = "lambda.1se",
    nfolds = 5,
    n_bootstrap = 100,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    min_bootstrap_success_fraction = 0.8,
    min_cv_rsq = 0,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    min_detected_cells_per_condition = 10,
    min_detection_fraction_per_condition = 0.01,
    seed = 12345L
  ),
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 30,
    rna_dims = 1:30,
    atac_dims = 2:30,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    seed = 12345L,
    overwrite = FALSE
  ),
  layer1_args = list(
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
  )
)
```

## Bootstrap provenance

Detailed Stage 1 outputs include:

```text
bootstrap_method
bootstrap_resampling_unit
bootstrap_sample_col
n_bootstrap_samples_total
min_bootstrap_samples_per_condition
bootstrap_fallback_reason
selection_frequency
sign_stability
```

These fields are propagated through `tf_peak_gene_condition_all`, `stability_diagnostics`, `celltype_fit_status`, and `group_status` before Stage 3 builds condition targets and complete-GPR cores.

## Console-only timing

Every public workflow stage prints elapsed time and final status in R:

```text
RegCompass timing: grn [success] 00:12:34.567
```

Timing is not stored in `result$timing`, stage-object `timing` fields, `step_timing.tsv`, or `00_execution_timing.tsv`.

## Compact result and detailed checkpoints

Primary compact tables include:

```r
result$table_manifest
result$reaction_ranking
result$condition_contrast
result$active_regulatory_edges
result$condition_target_genes
result$core_reactions
result$meta_module_summary
result$reaction_catalog
result$reaction_evidence
result$stage_provenance
```

Full candidates, coefficients, bootstrap diagnostics, metacell matrices, Layer 1 matrices, and reaction memberships remain in detailed stage checkpoints. Metacell-level tests are descriptive within-dataset analyses and are not biological-replicate treatment inference.

## Tutorials

1. [Quick start](docs/tutorial-01-quick-start.md)
2. [Stepwise audit](docs/tutorial-02-stepwise-audit.md)
3. [Restart and sensitivity](docs/tutorial-03-advanced-restart.md)
4. [Targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
5. [Condition differential analysis](docs/tutorial-05-condition-differential-analysis.md)
6. [Sample-aware bootstrap contract](docs/sample-aware-bootstrap.md)
