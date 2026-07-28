# Tutorial Level 2: stepwise run and Pando audit

Use this workflow when every stage should be saved and inspected independently.
The examples target RegCompassR 1.9.3 and Pando 1.2.1.

## Configure Stage 1 and Stage 4 workers

A supplied `BiocParallelParam` is forwarded to Pando. When `parallel = TRUE`
and `BPPARAM = NULL`, RegCompass instead enables Pando's native map backend.
When `parallel = FALSE`, Stage 1 is serial. Do not set
`pando_infer_args$parallel` separately.

```r
library(BiocParallel)

upstream_workers <- 6L
layer2_workers <- 30L

upstream_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(
    workers = upstream_workers,
    type = "SOCK",
    progressbar = TRUE
  )
} else {
  MulticoreParam(
    workers = upstream_workers,
    progressbar = TRUE
  )
}

layer2_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(
    workers = layer2_workers,
    type = "SOCK",
    progressbar = TRUE
  )
} else {
  MulticoreParam(
    workers = layer2_workers,
    progressbar = TRUE
  )
}
```

`BPPARAM = TRUE` is invalid because it is not a `BiocParallelParam` object.

## Stage 1: condition-comparable Pando inference

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 100L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      method = "shared_design_independent",
      candidate_screen = "motif_domain",
      tf_cor = 0.1,
      peak_cor = 0.01,
      alpha = 0.5,
      condition_mix = 1,
      condition_weight = "equal",
      reference_condition = "Control",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE
    )
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

Alternative without an explicit BiocParallel backend:

```r
step1_native <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn_native",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      reference_condition = "Control"
    )
  ),
  parallel = TRUE,
  BPPARAM = NULL
)
```

The resolved execution route is persisted:

```r
step1$params$pando_parallel
```

### Argument ownership

RegCompass forwards:

```text
pando_initiate_args → Pando::initiate_grn()
pando_motif_args    → Pando::find_motifs()
pando_infer_args    → Pando::infer_condition_grn()
```

The following are controlled by RegCompass and cannot be overridden inside the
nested lists:

```text
object, peak_assay, rna_assay, pfm, genome,
cell_type_col, condition_col, genes, network_name,
min_cells_per_condition, on_small_condition, BPPARAM
```

`aggregate_rna_col` and `aggregate_peaks_col` are rejected because Stage 1 fits
paired single cells and Stage 2 owns metacell aggregation.

### Candidate policy

The canonical default is `candidate_screen = "motif_domain"`. It retains the
structural motif/domain edge dictionary and lets elastic net select the
`TF RNA × peak ATAC` interaction predictors. `condition_union` and `pooled`
remain explicit marginal-screen sensitivity modes.

### Audit the Pando fit contract

```r
step1$grn_result$condition_fit_status
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition
step1$grn_result$tf_peak_gene_condition_effect

all_effects <- step1$grn_result$tf_peak_gene_condition_effect_all
all_effects[, c(
  "edge_id", "Group", "cell_type",
  "condition_estimate", "reference_estimate", "condition_effect",
  "eligible_in_condition", "comparable_to_reference"
)]

table(all_effects$comparable_to_reference, useNA = "ifany")
```

For edge `e`, condition `c`, and reference `r`:

```text
comparison_mask[e, c] = eligibility_mask[e, c] && eligibility_mask[e, r]
```

Only rows with `comparable_to_reference = TRUE` can enter the active
condition-effect table and downstream penalty projection.

The complete fit and transforms are written to:

```text
RegCompass_steps/01_grn/pando_condition_grn_fits.rds
RegCompass_steps/01_grn/pando_edge_predictor_transforms.tsv.gz
RegCompass_steps/01_grn/pando_tf_peak_gene_condition_effect_all.tsv.gz
RegCompass_steps/01_grn/pando_tf_peak_gene_condition_effect_active.tsv.gz
```

## Mouse Stage 1

Pando's bundled regulatory-region set is hg38. Mouse input therefore requires a
user-supplied build-matched `GRanges` object:

```r
library(BSgenome.Mmusculus.UCSC.mm10)

mouse_regions <- readRDS("mm10_regulatory_regions.rds")
stopifnot(methods::is(mouse_regions, "GenomicRanges"))

mouse_step1 <- rc_regcompass_step_grn(
  object = A_mouse,
  gem = mouse_gem,
  outdir = "RegCompass_mouse/01_grn",
  genome = BSgenome.Mmusculus.UCSC.mm10,
  species = "mouse",
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    pando_initiate_args = list(regions = mouse_regions),
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      reference_condition = "Control"
    )
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

The regulatory-region build must match the ATAC coordinates and motif-scanning
genome.

## Stage 2: condition-level multimodal metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  condition_col = "Group",
  celltype_col = "cell_type",
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  )
)

step2$pooled$metacell_meta
step2$pooled$cache_contract$analysis_args
```

Changing cells, assays, reductions, dimensions, embedding values, seed, gamma,
or thresholds invalidates the Stage 2 checkpoint. Set `overwrite = TRUE` to
rebuild.

## Stage 3: complete-GPR biological catalogue

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules"
)

step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$merged_modules$merged_core_reactions
```

Stage 3 performs one ordered expansion pass:

```text
complete-GPR cores
→ reactions in core-reaction subsystems
→ direct KEGG/Reactome equivalents
→ direct master-Rhea equivalents
```

It does not run FASTCORE.

## Stage 4: integrated reaction support

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp
)

step4$capacity_params$and_method
step4$evidence_formula
```

Stage 4 reconstructs Pando's standardized interaction predictor, applies the
stored condition-versus-reference coefficient, and excludes non-comparable
edges. It does not refit Pando or normalize coefficient effects by their
absolute sum.

## Stage 5: shared medium-specific model and LP scores

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_layer2",
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
  parallel = TRUE,
  BPPARAM = layer2_bp
)

step5$model_cache_summary
```

## Stage 6: assemble results

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_results",
  species = "human"
)

result$reaction_ranking
result$condition_contrast
result$grn$normalization_policy
```

See the complete [public function and API index](functions.md) for stage
arguments, return contracts, and restart boundaries.
