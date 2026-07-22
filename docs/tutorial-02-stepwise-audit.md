# Tutorial Level 2: stepwise run with audit gates

Use this level when every stage must be inspected before the next stage starts. It follows the same analysis design as Level 1 but exposes the intermediate objects and files.

Complete the installation and input checks in [Level 1](tutorial-01-quick-start.md) first. Use [Level 3](tutorial-03-advanced-restart.md) for restart, alternative media, solver selection, and resource controls.

## Stage map

| Stage | Primary input | Primary output | Required gate |
|---|---|---|---|
| 1. GRN | single-cell Seurat object, GEM, motifs, genome | condition × cell-type Pando GRNs | every required group is `ok` and has significant edges |
| 2. Metacells | original Seurat object | condition-only metacells with post hoc cell type | no ambiguous dominant-cell-type ties |
| 3. Meta-modules | Stages 1-2, GEM | core reactions and expanded modules | GRN/metacell coverage is complete; core reactions exist |
| 4. Layer 1 | metacells, modules, GEM | reaction-expression matrix | columns align exactly to metacell metadata |
| 5. Layer 2 | Layer 1, global module, medium | directional LP scores | targets were evaluated and feasible targets exist |
| 6. Results | Stages 1-5 | rankings and condition contrasts | outputs retain condition-specific and global modules |

## Common setup

```r
library(RegCompassR)
library(Pando)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

condition_col <- "dataset"
celltype_col <- "epithelial_or_stem"

gem <- rc_prepare_gem(species = "human", version = "2.0.0")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "high_glucose",
  species = "human"
)
rc_validate_gem(gem)
```

The same paired-cell Seurat object `A` is passed independently to Stages 1 and 2. Do not pass the internally normalized Stage 1 object into Stage 2.

## Stage 1: condition × cell-type single-cell GRNs

### Input

- original paired-cell RNA+ATAC object `A`;
- validated `gem`;
- `motifs`, not `motif2tf`;
- genome matching ATAC peak coordinates;
- complete condition and cell-type metadata.

### Run

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = condition_col,
  celltype_col = celltype_col,
  pando_args = list(
    min_cells = 100,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr"
    )
  )
)
```

RNA is normalized across all cells. ATAC TF-IDF uses all conditions within each cell type as the shared reference. Pando is then fitted independently for every condition × cell-type subset.

### Gate before Stage 2 or 3

```r
status1 <- step1$grn_result$sample_status
status1

stopifnot(
  all(status1$status == "ok"),
  all(status1$n_significant_edges > 0),
  nrow(step1$grn_result$tf_peak_gene_significant) > 0
)

step1$grn_result$pando_installation
head(step1$grn_result$tf_peak_gene_significant)
```

Do not continue when a required group failed, was skipped for too few cells, or has no significant edge.

### Files

- `single_cell_grn.rds`;
- `step_grn.rds`;
- `pando_group_status.tsv.gz`;
- `pando_tf_peak_gene_all.tsv.gz`;
- `pando_tf_peak_gene_significant.tsv.gz`;
- optional `pando_objects/*.rds`.

## Stage 2: condition-only metacells

### Input

The original `A`. Condition is the only varying metacell stratum. Cell type and sample are not metacell grouping variables.

### Run

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  condition_col = condition_col,
  celltype_col = celltype_col,
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 75,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  )
)
```

Each metacell receives a dominant member-cell type after construction. Purity, mixed-cell-type status, and the full composition are retained. An exact dominant-cell-type tie stops the workflow.

### Gate before Stage 3

```r
meta2 <- step2$pooled$metacell_meta

stopifnot(
  nrow(meta2) > 0,
  !anyNA(meta2[[condition_col]]),
  !anyNA(meta2[[celltype_col]]),
  !any(meta2$dominant_celltype_tied %in% TRUE),
  setequal(colnames(step2$metacell_object), meta2$metacell_id)
)

table(meta2[[condition_col]], meta2[[celltype_col]])
summary(meta2$dominant_celltype_fraction)
```

### Files

- `step_metacells.rds`;
- `merged_metacell_object.rds`;
- `metacell_metadata.tsv.gz`;
- `metacell_membership.tsv.gz`;
- `metacell_celltype_composition.tsv.gz`;
- `metacell_celltype_summary.tsv.gz`.

## Stage 3: core reactions and meta-modules

### Input

Stage 1, Stage 2, and the same GEM used for Stage 1.

### Run

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules",
  layer1_args = list(
    top_k_neighbors = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    local_fastcore = TRUE,
    local_fastcore_args = list(
      solver = "highs",
      strict = TRUE
    )
  )
)
```

The stage:

1. validates bidirectional condition × cell-type coverage between GRNs and metacells;
2. projects the Pando network onto metabolic genes;
3. maps complete-GPR core reactions;
4. adds reactions from core-reaction subsystems;
5. adds reactions sharing KEGG, Reactome, or identical master-Rhea identifiers;
6. adds only the FASTCORE support reactions needed for local feasibility.

### Gate before Stage 4

```r
coverage3 <- step3$group_coverage
core3 <- step3$condition_modules$core_gene_reaction
membership3 <- step3$condition_modules$reaction_membership

stopifnot(
  all(coverage3$coverage_complete),
  nrow(core3) > 0,
  nrow(membership3) > 0,
  all(core3$reaction_id %in% colnames(gem$S)),
  all(membership3$reaction_id %in% colnames(gem$S))
)

coverage3
head(core3)
head(membership3)
step3$condition_modules$local_fastcore_summary
```

### Files

- `grn_metacell_group_coverage.tsv.gz`;
- `core_gene_reaction.tsv.gz`;
- `meta_module_reactions.tsv.gz`;
- `condition_meta_modules.rds`;
- `global_meta_modules.rds`;
- `local_fastcore/` diagnostics and optional models.

## Stage 4: RNA+ATAC reaction expression

### Run

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  tau = 0.20
)
```

### Gate before Stage 5

```r
stopifnot(
  nrow(step4$reaction_expression) > 0,
  ncol(step4$reaction_expression) > 0,
  identical(
    colnames(step4$reaction_expression),
    as.character(step4$unit_meta$pool_id)
  ),
  all(is.finite(step4$reaction_expression))
)

dim(step4$reaction_expression)
head(step4$gpr_diagnostics)
```

`step_layer1.rds` contains RNA support, the ATAC regulatory modifier, integrated gene support, GPR diagnostics, and reaction expression.

## Stage 5: directional COMPASS-like scoring

### Run

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
    time_limit = 60
  )
)
```

The selected solver is checked before medium-constrained model construction. A missing solver package is reported separately from true biological infeasibility.

### Gate before Stage 6

```r
stopifnot(
  any(step5$evaluated),
  nrow(step5$lp_diagnostics) > 0,
  nrow(step5$model_cache_summary) > 0
)

table(step5$evaluated)
table(step5$feasible)
head(step5$lp_diagnostics)
step5$model_cache_summary
```

A run can contain biologically blocked target directions, but a result with no evaluated target or no feasible target requires inspection before interpretation.

### Files

- `step_layer2.rds`;
- exported score, penalty, `vmax`, feasibility, and diagnostic tables;
- persistent `model_cache/`.

## Stage 6: final result assembly

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_results"
)
```

### Final gate

```r
stopifnot(
  nrow(result$reaction_ranking) > 0,
  identical(result$version, "1.8.1"),
  !is.null(result$condition_grn_meta_modules),
  !is.null(result$global_grn_meta_modules)
)

head(result$reaction_ranking)
head(result$condition_summary)
head(result$condition_contrast)
```

Final files are `step_comparison.rds` and `regcompass_result.rds`.

## Stop conditions

Stop and diagnose rather than continuing when any of these occur:

- a condition × cell-type Pando group is not `ok`;
- a metacell has an ambiguous dominant-cell-type tie;
- GRN/metacell group coverage is incomplete;
- no complete-GPR core reaction remains;
- Layer 1 reaction-expression columns do not align to metacell metadata;
- Layer 2 evaluated no target;
- solver installation is reported as missing.
