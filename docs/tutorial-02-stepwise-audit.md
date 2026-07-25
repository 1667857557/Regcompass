# Tutorial Level 2: stepwise run

Use this tutorial when each RegCompass stage should be run and saved independently.

## Stage 1: infer condition-by-cell-type Pando evidence

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  pfm = motif2tf,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 100,
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    )
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp
)

step1$grn_result$sample_status
step1$grn_result$target_metabolic_genes
```

The candidate targets are all Human-GEM GPR genes present in the RNA assay. Unless `pando_args$pando_initiate_args$regions` is supplied, the human workflow uses:

```r
data("phastConsElements20Mammals.UCSC.hg38", package = "Pando")
data("SCREEN.ccRE.UCSC.hg38", package = "Pando")
regions <- union(
  phastConsElements20Mammals.UCSC.hg38,
  SCREEN.ccRE.UCSC.hg38
)
```

This union is passed to `Pando::initiate_grn(regions = regions)`. The default is hg38-specific; mouse or other genomes require an explicit region object.

## Stage 2: construct condition-level metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  sample_col = NULL,
  condition_col = "Group",
  celltype_col = "cell_type",
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 30,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  )
)

step2$pooled$metacell_meta
```

## Stage 3: construct complete-GPR biological meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules",
  meta_module_args = list(
    expansion_mode = "ordered_once"
  )
)
```

Stage 3 no longer projects targets through shared TFs and does not calculate GRN connected components. For each `condition × cell type`, it performs the following operations:

```text
significant Pando TF–peak–Human-GEM-target rows
→ unique supported metabolic target genes
→ complete-GPR core reactions
→ core-subsystem expansion
→ KEGG/Reactome reaction-equivalence expansion
→ master-Rhea reaction-equivalence expansion
→ biological meta-module
```

A positive or negative Pando coefficient both count as regulatory evidence. A reaction is core only when at least one complete GPR AND branch is contained in the supported target-gene set. Partial complexes remain diagnostic and do not anchor expansion.

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
table(step3$condition_modules$reaction_membership$inclusion_stage)

catalogue <- step3$merged_modules
catalogue$merged_core_reactions
catalogue$merged_reaction_membership
```

## Stage 4: calculate integrated RNA+ATAC reaction support

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

`gpr_and_method` accepts the COMPASS functions `"min"`, `"median"`, and `"mean"`. RegCompass defaults to `"min"`, so the least-supported required subunit limits a multi-gene GPR complex. The canonical isozyme OR rule remains additive. The former Boltzmann soft-min and `tau` API have been deleted.

The selected rule is recorded in:

```r
step4$capacity_params$and_method
step4$evidence_formula
```

## Stage 5: build the medium-constrained model and score reactions

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

Stage 5 first applies the selected medium and performs global FASTCORE completion to construct the union GEM. `completion_time_limit` applies only to this construction phase. The completed union GEM is then cached and reused for directional scoring; scoring LPs have no time-limit parameter. See [medium presets](medium-presets.md) for available presets and custom media.

```r
step5$model_cache_summary[, c(
  "medium_scenario",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "n_reactions",
  "file"
)]
```

## Stage 6: assemble annotated results

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
result$condition_grn_meta_modules$supported_metabolic_genes
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
```
