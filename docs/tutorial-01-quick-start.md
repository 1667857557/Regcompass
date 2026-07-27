# Tutorial Level 1: minimal one-shot run

**This tutorial:** runs the complete RegCompassR 1.8.10 workflow, uses donor/sample-aware Stage 1 bootstrap when available, and introduces the compact final result.

**Next:** [Tutorial 2 — stepwise run and audit](tutorial-02-stepwise-audit.md) reproduces the same analysis stage by stage and exposes the detailed objects intentionally omitted from `result`.

## Workflow

```text
shared GREAT-domain Pando TF–peak–target structural candidates per cell type
→ condition-aware TF×peak/target observability filter
→ one shared model edge universe across conditions
→ condition-balanced direct condition-theta elastic net
→ derived global backbone + zero-sum condition deviations
→ condition-stratified sample/donor cluster bootstrap when sample_col is valid
→ explicit warning and condition-stratified cell-bootstrap fallback otherwise
→ condition-specific bootstrap-active sub-GRNs
→ condition target genes
→ complete-GPR reaction cores
→ biological module expansion
→ shared medium-specific union GEM
→ RNA+ATAC penalties
→ direction-specific LP scoring
→ compact analysis result
```

## 1. Prepare the model

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

The paired Seurat object must contain RNA and ATAC assays, the requested PCA/LSI reductions, and complete metadata columns named by:

```text
condition_col
celltype_col
```

A biological sample/donor column is optional but recommended:

```text
sample_col
```

A valid `sample_col` is used only for Stage 1 bootstrap stability. Sample IDs are sampled with replacement separately within each condition, and every selected sample contributes all of its cells. When `sample_col = NULL`, or the named column does not exist, RegCompass prints the exact fallback reason and uses condition-stratified cell resampling. An existing sample column with missing or empty IDs is an error.

Stage 2 uses the current SuperCell2 contract:

```text
RegCompass splits cells by condition
→ SCimplify_for_Seurat(label = celltype_col)
→ exact label-preserving metacells
```

The sample column does not alter Stage 2 strata and does not create sample-specific metabolic models.

## 2. Run the complete workflow

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
    min_cells_per_stratum = 500,
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
  layer2_workers = 30,
  progress = TRUE
)
```

The canonical peak-to-gene method is `GREAT`. Its basal-plus-extension domains permit distal structural candidates without using target-expression correlation to admit edges. With `extend = 1000000`, the hypothesis space is deliberately broad and is subsequently constrained by condition-aware observability, regularisation, CV, and bootstrap stability.

The canonical default is `n_bootstrap = 100L`. RegCompass reports the Monte Carlo standard error and Wilson 95% interval of each condition-edge selection frequency.

With a valid sample column, the Stage 1 method is:

```text
condition_stratified_sample_cluster_nonparametric
```

If the sample column is omitted or absent, R prints the reason and the method becomes:

```text
condition_stratified_cell_nonparametric_fallback
```

Each public stage also prints elapsed time and final status after its final artifact is committed, for example:

```text
RegCompass timing: grn [success] 00:12:34.567
```

Timing is not stored in `result$timing`, stage-object `timing` fields, `step_timing.tsv`, or `00_execution_timing.tsv`.

## 3. Understand the key calculations

For candidate edge `e = (TF t, peak p, target g)`:

\[
x_{e,u}=T_{t,u}A_{p,u}.
\]

The pooled Pando detection thresholds remain zero so a condition-restricted regulator is not removed before multitask fitting. RegCompass retains an edge in the shared model universe only when the non-zero `TF × peak` predictor and target RNA each occur in at least

\[
m_c=\min\left(n_c,\max\left(10,\left\lceil0.01n_c\right\rceil\right)\right)
\]

cells of at least one condition.

RegCompass directly estimates

\[
y_u^\circ=\sum_e\widetilde x_{e,u}\theta_{e,c(u)}+\varepsilon_u.
\]

The L1 penalty acts on each \(\theta_{e,c}\), allowing an edge to be exactly zero in one condition and non-zero in another. The global backbone and deviations are derived summaries:

\[
\beta_e=\frac1C\sum_c\theta_{e,c},
\qquad
\delta_{e,c}=\theta_{e,c}-\beta_e,
\qquad
\sum_c\delta_{e,c}=0.
\]

`alpha = 0.5` combines sparse selection with ridge stabilization of correlated TF–peak predictors. `global_penalty_factor` and `deviation_penalty_factor` are compatibility aliases for one common direct-theta penalty and must be equal.

When a valid sample column exists, each bootstrap resamples the observed number of sample/donor clusters with replacement within every condition and includes all cells from each selected cluster. Cluster sizes are preserved, so bootstrap cell counts can vary. The fallback instead resamples each condition at its original cell count. Both modes recalculate within-condition centring and fit at the full-data selected lambda.

An edge becomes active only after passing selection-frequency, sign-stability, effect-size, strictly positive out-of-fold CV R-squared, and bootstrap-completion thresholds.

For sign stability

\[
\rho=|2q-1|,
\]

so `min_sign_stability = 0.8` requires at least 90% agreement on one sign among selected bootstrap fits.

For condition target-gene set `G_c`, a reaction is core only when one complete GPR branch is present:

\[
Core_{r,c}=1\iff\exists k:B_{r,k}\subseteq G_c.
\]

Positive and negative active edges both establish regulated-gene membership; the coefficient sign is retained in the ATAC projection.

See [Pando and multitask GRN parameter policy](grn-parameter-policy.md) and [Sample-aware bootstrap contract](sample-aware-bootstrap.md) for the full default rationale, resampling semantics, and supported sensitivity analyses.

## 4. Inspect the compact final result

```r
result$schema_version
result$version
result$table_manifest
```

Primary tables are:

```r
result$reaction_ranking
result$condition_contrast
result$active_regulatory_edges
result$condition_target_genes
result$core_reactions
result$meta_module_summary
result$grn_metacell_group_coverage
result$reaction_catalog
result$reaction_evidence
```

`reaction_ranking` and `condition_contrast` contain only analysis columns. Reaction names, formulas, GPRs, and database cross-references occur once in `reaction_catalog` rather than being repeated in every ranking row.

`active_regulatory_edges` retains the bootstrap method, resampling unit, sample column, sample counts, and fallback reason. The global Stage 1 policy is also retained in:

```r
result$params$sample_col
result$params$bootstrap_resampling_unit
result$params$bootstrap_fallback_reason
result$stage_provenance$bootstrap_policy
```

The final object retains `result$microcompass` because direction-specific condition testing and plotting require the unit-level scores. It does not embed the full GRN candidate universe, all coefficient rows, metacell assay matrices, Layer 1 matrices, or full module membership.

```r
result$stage_provenance$detailed_sources
```

## 5. Locate detailed stage outputs

The one-shot run writes stage checkpoints under the output directory. Use them only when detailed auditing is required:

```r
step1 <- readRDS("RegCompass_result/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_result/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_result/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_result/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_result/05_layer2/step_layer2.rds")
```

Actual subdirectory names may follow the one-shot stage layout shown in the run log; `result$stage_provenance` records which scientific information belongs to each stage.

Audit the Stage 1 sampling provenance directly:

```r
step1$grn_result$bootstrap_policy
head(step1$grn_result$stability_diagnostics[, intersect(c(
  "bootstrap_method",
  "bootstrap_resampling_unit",
  "bootstrap_sample_col",
  "n_bootstrap_samples_total",
  "min_bootstrap_samples_per_condition",
  "bootstrap_fallback_reason",
  "selection_frequency",
  "sign_stability"
), colnames(step1$grn_result$stability_diagnostics)), drop = FALSE])
```

These fields are finalized before Stage 3 consumes `tf_peak_gene_significant`, so condition target genes and complete-GPR cores are derived from the correctly resampled stability-selected sub-GRN.

## 6. Exit checks before Tutorial 2 or Tutorial 5

```r
stopifnot(
  identical(result$version, "1.8.10"),
  identical(result$params$sample_col, "sample_id"),
  identical(result$params$bootstrap_resampling_unit, "sample"),
  isTRUE(!result$stage_provenance$detailed_intermediates_embedded),
  nrow(result$reaction_ranking) > 0L,
  nrow(result$reaction_catalog) > 0L
)
```

For a transparent stage-by-stage reconstruction, continue to [Tutorial 2](tutorial-02-stepwise-audit.md). For biological interpretation after a completed run, continue to [Tutorial 5](tutorial-05-condition-differential-analysis.md).
