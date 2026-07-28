# RegCompassR

RegCompassR 1.8.10 implements a shared-background regulatory–metabolic workflow for paired single-cell RNA+ATAC data.

## Canonical architecture

```text
all conditions within one cell type
→ one validated GREAT-domain Pando structural TF–peak–metabolic-gene universe
→ condition-aware TF×peak/target observability filter
→ one shared model edge universe across conditions
→ condition-balanced direct condition-theta elastic net
→ derived cross-condition GRN backbone + zero-sum deviations
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

The default `grn_mode = "multitask_shared_backbone"` fits condition-specific coefficients directly. For edge \(e=(TF,peak,target)\),

\[
\beta_e=\frac1C\sum_c\theta_{e,c},
\qquad
\delta_{e,c}=\theta_{e,c}-\beta_e,
\qquad
\sum_c\delta_{e,c}=0.
\]

All conditions of one cell type use the same structural edge dictionary, filtered model edge universe, predictor scale, penalty structure, and candidate ordering. The elastic-net L1 penalty acts directly on \(\theta_{e,c}\), so an edge may be exactly zero in one condition and non-zero in another. A reaction becomes a condition core only when at least one complete GPR branch is contained in the condition target-gene set.

## Installation

The validated default profile remains Seurat v4:

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

Pando 1.1.3 or later is required. RegCompass validates the Pando version-2 design fingerprint before fitting the multitask GRN. SeuratObject/Seurat 5.x with Signac 1.12–1.x is also accepted; see [Seurat compatibility](docs/seurat-compatibility.md).

## Metadata and bootstrap contract

The canonical workflow requires `condition_col` and `celltype_col`. A biological donor or sample column is optional but strongly recommended through `sample_col`.

`sample_col` is accepted by `rc_run_regcompass()`, `rc_run_regcompass_one_shot()`, and `rc_regcompass_step_grn()` only. With a valid column, RegCompass performs condition-stratified sample-cluster bootstrap: sample IDs are drawn with replacement separately inside each condition, and every selected donor/sample contributes all of its cells.

When `sample_col = NULL`, or the named metadata column does not exist, RegCompass prints a warning containing the exact fallback reason and performs condition-stratified cell resampling instead. An existing sample column containing missing or empty identifiers is rejected because silently assigning such cells to pseudo-donors would be invalid.

The Stage 1 outputs retain:

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

Sample-cluster bootstrap estimates reproducibility across the observed biological samples. Cell fallback estimates cell-resampling reproducibility only. Neither converts metacells into biological-replicate inference or provides formal treatment-effect inference.

See [Sample-aware bootstrap contract](docs/sample-aware-bootstrap.md).

## SuperCell2 contract

Stage 2 remains condition based and is deliberately independent of the bootstrap sampling unit:

```text
RegCompass splits the Seurat object by strata_cols = condition_col
→ calls SCimplify_for_Seurat(label = celltype_col)
→ validates exact label-pure membership
```

There is no artificial condition-pool metadata field. In the Stage 2 public API, `sample_col` is intentionally absent: donor/sample labels control Stage 1 resampling but do not create sample-specific metacells or sample-specific metabolic models.

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
  ),
  upstream_workers = 6,
  layer2_workers = 30
)
```

The canonical peak-to-gene rule is `GREAT`. With `extend = 1000000`, it defines a broad distal regulatory-domain hypothesis without using target-expression correlation to admit candidates. `Signac` remains available as an explicit sensitivity analysis.

The canonical default is `n_bootstrap = 100L`. RegCompass reports bootstrap `selection_frequency`, `sign_stability`, Monte Carlo standard errors, and Wilson 95% intervals. Finite `max_edges_per_target` values are rejected because Pando candidate order is deterministic but is not an evidence ranking.

## Stepwise workflow and output chaining

The inspectable stages are:

1. `rc_regcompass_step_grn()` — shared structural candidates, multitask coefficients, CV, bootstrap stability, and sampling provenance.
2. `rc_regcompass_step_metacells()` — condition-stratified, cell-type-labelled SuperCell2 metacells.
3. `rc_regcompass_step_meta_modules()` — bootstrap-active targets, complete-GPR cores, and biological reaction modules.
4. `rc_regcompass_step_layer1()` — RNA support, ATAC modifier, GPR aggregation, and reaction penalties.
5. `rc_regcompass_step_layer2()` — shared medium-specific union GEM and directional LP scoring.
6. `rc_regcompass_step_results()` — compact result tables and provenance.
7. `rc_regcompass_step_target_union()` — optional direct-equivalent second-pass target scoring in the cached model.

Each stage validates the class, workflow parameters, condition-by-cell-type coverage, GEM fingerprint, reaction identifiers, and matrix/unit alignment received from the previous stage. Detailed intermediate matrices and candidate tables remain in the detailed stage checkpoints rather than being duplicated in the final object.

The compact final result contains:

```r
result$table_manifest
result$reaction_ranking
result$condition_contrast
result$active_regulatory_edges
result$condition_target_genes
result$core_reactions
result$meta_module_summary
result$grn_metacell_group_coverage
result$reaction_catalog
result$reaction_evidence
result$stage_provenance
```

`active_regulatory_edges` retains bootstrap method, resampling unit, sample column, sample counts, and fallback reason. The global Stage 1 bootstrap policy is retained in both `result$params` and `result$stage_provenance`. Reaction annotations are stored once in `reaction_catalog`; condition- and cell-type-specific RNA versus RNA+ATAC provenance is stored in `reaction_evidence`.

## Execution timing

Every public workflow stage prints elapsed time and final status to the R console after its final artifact has been committed:

```text
RegCompass timing: grn [success] 00:12:34.567
```

Timing is not stored in `result$timing`, stage-object `timing` fields, `step_timing.tsv`, or `00_execution_timing.tsv`. One-shot execution also removes stale timing files from an existing output directory.

## Compare the same reaction across conditions

Use `rc_test_condition_reactions()` to compare a fixed reaction, direction, medium, and cell type across conditions. With three or more conditions, RegCompass returns a Kruskal-Wallis omnibus test and pairwise Wilcoxon tests.

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  cell_types = "stem-cell_like",
  min_units = 5,
  p_adjust_method = "BH",
  p_adjust_scope = "celltype_contrast_medium",
  outdir = "RegCompass_result/07_condition_statistics"
)

condition_stats$omnibus
condition_stats$pairwise
```

Plot one direction of one reaction:

```r
p <- rc_plot_condition_reaction(
  result,
  reaction_id = "MAR06231",
  cell_type = "stem-cell_like",
  target_direction = "reverse",
  medium_scenario = "high_glucose",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  annotation_p = "p_adj"
)
print(p)
```

The plot shows one point per metacell, overlays the condition distributions, and adds significance brackets for requested contrasts. These are descriptive metacell-level comparisons; metacells are not independent biological replicates. Interpret adjusted P values together with effect sizes, GPR structure, and `reaction_evidence`.

See [Condition-associated reaction statistics](docs/condition-reaction-statistics.md) for reaction formulas, evidence classes, gene-based selection, multi-condition tests, and plotting details.

## Legacy mode

Independent condition-by-cell-type Pando fitting remains available for reproducibility:

```r
result_legacy <- rc_run_regcompass_one_shot(
  ...,
  grn_mode = "legacy_condition_pando",
  pando_args = list(
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    ),
    padj_threshold = 0.05,
    min_model_rsq = 0.1
  )
)
```

Legacy adjusted-P-value thresholds are not interpreted as significance thresholds in multitask mode.

## Tutorials

1. [Quick start](docs/tutorial-01-quick-start.md)
2. [Stepwise audit](docs/tutorial-02-stepwise-audit.md)
3. [Restart and sensitivity](docs/tutorial-03-advanced-restart.md)
4. [Targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
5. [Condition analysis from GRN to reaction](docs/tutorial-05-condition-differential-analysis.md)

Additional references:

- [Condition-associated reaction statistics](docs/condition-reaction-statistics.md)
- [Sample-aware bootstrap contract](docs/sample-aware-bootstrap.md)
- [GRN parameter policy](docs/grn-parameter-policy.md)
- [Shared-background multitask GRN mathematics](docs/multitask-shared-grn.md)
- [Stage contracts](docs/stage-interface-contracts.md)
- [Medium presets](docs/medium-presets.md)
- [Direction-aware reporting](docs/direction-aware-condition-reporting.md)
- [Seurat compatibility](docs/seurat-compatibility.md)
