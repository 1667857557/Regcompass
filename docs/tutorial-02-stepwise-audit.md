# Tutorial Level 2: stepwise run and audit

Run each RegCompassR 1.8.8 stage independently when intermediate objects, tables, or restart points must be inspected.

## Parallel backends

```r
library(BiocParallel)

upstream_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 6L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 6L, progressbar = TRUE)
}

layer2_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 30L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 30L, progressbar = TRUE)
}
```

Stage 1 parallelises cell types and target models; each individual `glmnet` fit remains single-threaded. Stage 4 reuses `upstream_bp`; Stage 5 uses `layer2_bp`.

## Stage 1: shared GRN background and condition sub-GRNs

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  sample_col = "sample_id",
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
    n_stability = 50,
    stability_fraction = 0.8,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    seed = 12345L
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

Inspect the Stage 1 contract:

```r
step1$grn_result$celltype_fit_status
step1$grn_result$sample_status
step1$grn_result$tf_peak_gene_candidates
step1$grn_result$tf_peak_gene_global
step1$grn_result$tf_peak_gene_condition_all
step1$grn_result$tf_peak_gene_significant
step1$grn_result$condition_target_genes
step1$grn_result$target_model_diagnostics
step1$grn_result$stability_diagnostics
```

For a given cell type, every condition must have the same `edge_universe_id`. The condition coefficient fields satisfy:

```text
effective_estimate = global_estimate + condition_deviation
estimate = effective_estimate × selection_frequency × sign_stability
```

The sum of `condition_deviation` over all conditions is zero for each edge. `padj` is `NA` in multitask mode because regularised stability selection, rather than a classical coefficient test, defines active edges.

Stage 1 writes:

```text
pando_celltype_fit_status.tsv.gz
pando_group_status.tsv.gz
pando_tf_peak_gene_candidates.tsv.gz
pando_tf_peak_gene_global.tsv.gz
pando_tf_peak_gene_condition_all.tsv.gz
pando_tf_peak_gene_all.tsv.gz
pando_tf_peak_gene_significant.tsv.gz
condition_target_genes.tsv.gz
pando_target_model_diagnostics.tsv.gz
pando_edge_stability.tsv.gz
single_cell_grn.rds
step_grn.rds
```

## Stage 2: construct condition-level metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  sample_col = "sample_id",
  condition_col = "Group",
  celltype_col = "cell_type",
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
  )
)

step2$pooled$metacell_meta
step2$pooled$cache_contract$analysis_args
```

Metacells are built within condition, with cell type used as a SuperCell2 label and audited after aggregation. The sample column is retained as provenance; it does not split the Stage 2 strata.

## Stage 3: condition-specific complete-GPR cores

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules"
)
```

For each `condition × cell type`:

```text
stable active TF–peak–metabolic-gene edges
→ unique regulated metabolic target genes
→ complete-GPR core reactions
→ core-reaction subsystems
→ direct KEGG/Reactome equivalents
→ direct master-Rhea equivalents
→ biological condition meta-module
```

A reaction is core only when at least one complete GPR AND branch is present. Partial complexes remain diagnostic and cannot anchor expansion.

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
step3$merged_modules$source_edge_universe_ids
```

The merged output is a reaction catalogue, not a GEM.

## Stage 4: RNA+ATAC reaction support

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
```

Inspect:

```r
step4$gene_support_rna
step4$gene_regulatory_modifier
step4$gene_support_multiome
step4$reaction_expression
step4$gpr_diagnostics
```

The regulatory modifier uses only metacell ATAC during projection. TFs sharing the same measured peak are signed-summed, and one target denominator is shared across conditions. A gene without an active condition edge has modifier zero and exact RNA-only support.

`gpr_and_method` accepts `"min"`, `"median"`, or `"mean"`; the canonical default is `"min"`. Isozyme OR branches remain additive.

## Stage 5: shared medium-specific union GEM and LP scoring

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
```

For each medium, Stage 5 performs the only FASTCORE completion and saves one union GEM. All conditions and metacells reuse the exact same reaction IDs, stoichiometric matrix and bounds.

```r
step5$model_cache_summary
step5$source_core_reactions
step5$source_merged_reaction_membership
step5$union_gem_policy
```

`completion_time_limit` limits union-GEM construction only; scoring LPs do not accept a time-limit parameter.

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
```

```r
result$schema_version
result$version
result$grn_mode
result$reaction_ranking
result$condition_contrast
result$reaction_catalog
result$reaction_evidence
```

Each public stage writes a `step_*.rds` restart object and validates the GEM fingerprint, metadata/assay signature, metacell order and upstream class before continuing.

See [multitask GRN mathematics and object contracts](multitask-shared-grn.md).
