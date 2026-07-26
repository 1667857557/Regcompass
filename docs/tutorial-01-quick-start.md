# Tutorial Level 1: minimal one-shot run

**This tutorial:** runs the complete RegCompassR 1.8.9 workflow and introduces the compact final result.

**Next:** [Tutorial 2 — stepwise run and audit](tutorial-02-stepwise-audit.md) reproduces the same analysis stage by stage and exposes the detailed objects intentionally omitted from `result`.

## Workflow

```text
shared Pando TF–peak–target candidates per cell type
→ condition-balanced multitask elastic net
→ global backbone + condition deviations
→ condition-stratified bootstrap-active sub-GRNs
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

No sample column is accepted or interpreted by the canonical workflow.

Stage 2 uses the current SuperCell2 contract:

```text
RegCompass splits cells by condition
→ SCimplify_for_Seurat(label = celltype_col)
→ exact label-preserving metacells
```

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
    min_bootstrap_success_fraction = 0.8,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
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

The package default is `n_bootstrap = 50L`; `100` is suitable for a final run with lower Monte Carlo error in empirical selection frequencies.

## 3. Understand the key calculations

For candidate edge `e = (TF t, peak p, target g)`:

\[
x_{e,u}=T_{t,u}A_{p,u}.
\]

The condition coefficient is

\[
\theta_{e,c}=\beta_e+\delta_{e,c},
\qquad \sum_c\delta_{e,c}=0.
\]

Each bootstrap resamples every condition with replacement at its original cell count and recalculates within-condition centring. An edge becomes active only after passing selection-frequency, sign-stability, effect-size, CV-reliability, and bootstrap-completion thresholds.

For condition target-gene set `G_c`, a reaction is core only when one complete GPR branch is present:

\[
Core_{r,c}=1\iff\exists k:B_{r,k}\subseteq G_c.
\]

Positive and negative active edges both establish regulated-gene membership; the coefficient sign is retained in the ATAC projection.

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

## 6. Exit checks before Tutorial 2 or Tutorial 5

```r
stopifnot(
  identical(result$version, "1.8.9"),
  isTRUE(!result$stage_provenance$detailed_intermediates_embedded),
  nrow(result$reaction_ranking) > 0L,
  nrow(result$reaction_catalog) > 0L
)
```

For a transparent stage-by-stage reconstruction, continue to [Tutorial 2](tutorial-02-stepwise-audit.md). For biological interpretation after a completed run, continue to [Tutorial 5](tutorial-05-condition-differential-analysis.md).
