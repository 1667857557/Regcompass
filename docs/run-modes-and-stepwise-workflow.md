# RegCompassR stepwise workflow

This is the canonical 1.8.2 workflow. Each section states the stage input, output, and minimum inspection required before continuing.

## Installation

```r
install.packages(c("remotes", "highs"))
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0", upgrade = "never")
remotes::install_github("1667857557/Pando_regcompass", upgrade = "never")
remotes::install_github("1667857557/Regcompass", upgrade = "never")
```

A locally downloaded Pando source archive is supported:

```r
install.packages("~/Pando_regcompass.tar.gz", repos = NULL, type = "source")
```

RegCompass checks the required Pando API. GitHub remote metadata are not required for a local installation.

## Input contract

The same paired-cell Seurat multiome object is passed to Steps 1 and 2.

| Input | Requirement |
|---|---|
| RNA | Raw-count assay with GEM-compatible gene symbols |
| ATAC | Peak-count `ChromatinAssay` |
| Cells | RNA, ATAC, and metadata contain the same cell IDs; order may differ |
| Metadata | Complete condition and cell-type columns |
| Genome | Matches the ATAC peak coordinates |
| `pfm` | Pando/motifmatchr-compatible PFM/PWM collection |
| Fragments | Not needed when `fragment_files = FALSE` |

Use `Pando::motifs` as `pfm`; do not pass the `motif2tf` annotation table.

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

condition_col <- "dataset"
celltype_col <- "epithelial_or_stem"

stopifnot(
  inherits(A, "Seurat"),
  all(c("RNA", "ATAC") %in% names(A@assays)),
  inherits(A[["ATAC"]], "ChromatinAssay"),
  all(c(condition_col, celltype_col) %in% colnames(A@meta.data)),
  !anyNA(A@meta.data[[condition_col]]),
  !anyNA(A@meta.data[[celltype_col]]),
  all(nzchar(trimws(as.character(A@meta.data[[condition_col]])))),
  all(nzchar(trimws(as.character(A@meta.data[[celltype_col]]))))
)

gem <- rc_prepare_gem(species = "human", version = "2.0.0")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "high_glucose",
  species = "human"
)
rc_validate_gem(gem)
```

RNA is normalized globally. ATAC TF-IDF is calculated once within each cell type across all conditions. A peak absent from one cell type remains zero and is not passed to that cell type's `RunTFIDF` call.

## Step 1: single-cell GRNs

**Input:** `A`, `gem`, `motifs`, genome, condition column, cell-type column.

**Output used by Step 3:** `step1$grn_result`.

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

step1$grn_result$sample_status
head(step1$grn_result$tf_peak_gene_significant)
step1$grn_result$pando_installation
```

Do not continue if any Pando group failed or if required condition × cell-type groups have no significant edges.

Files:

- `single_cell_grn.rds`, `step_grn.rds`
- `pando_group_status.tsv.gz`
- `pando_tf_peak_gene_all.tsv.gz`
- `pando_tf_peak_gene_significant.tsv.gz`
- optional `pando_objects/*.rds`

## Step 2: condition-only metacells

**Input:** original `A`, not the normalized object internal to Step 1.

**Output used by Steps 3-4:** `step2$pooled` and `step2$metacell_object`.

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

table(
  step2$pooled$metacell_meta[[condition_col]],
  step2$pooled$metacell_meta[[celltype_col]]
)
summary(step2$pooled$metacell_meta$dominant_celltype_fraction)
```

SuperCell2 is grouped only by condition. Cell type is assigned afterwards from member cells. Exact dominant-cell-type ties stop the workflow.

Files:

- `step_metacells.rds`, `merged_metacell_object.rds`
- `metacell_metadata.tsv.gz`
- `metacell_membership.tsv.gz`
- `metacell_celltype_composition.tsv.gz`
- `metacell_celltype_summary.tsv.gz`

## Step 3: core reactions and meta-modules

**Input:** Step 1, Step 2, and the same GEM.

**Output used by Steps 4-5:** `step3$condition_modules` and `step3$global_modules`.

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

step3$group_coverage
head(step3$condition_modules$core_gene_reaction)
head(step3$condition_modules$reaction_membership)
step3$condition_modules$local_fastcore_summary
```

The stage validates GRN ↔ metacell group coverage, maps complete-GPR core reactions, expands through core subsystems and shared KEGG/Reactome/master-Rhea identifiers, then adds local FASTCORE support reactions.

Files:

- `grn_metacell_group_coverage.tsv.gz`
- `core_gene_reaction.tsv.gz`
- `meta_module_reactions.tsv.gz`
- `condition_meta_modules.rds`, `global_meta_modules.rds`
- `local_fastcore/` diagnostics and optional models

## Step 4: reaction expression

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  tau = 0.20
)

stopifnot(identical(
  colnames(step4$reaction_expression),
  step4$unit_meta$pool_id
))
dim(step4$reaction_expression)
```

`step_layer1.rds` contains RNA support, ATAC modifiers, integrated gene support, GPR diagnostics, and reaction expression.

## Step 5: directional scoring

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

table(step5$feasible)
head(step5$lp_diagnostics)
step5$model_cache_summary
```

The selected solver is checked before model construction. A missing solver package is reported as an installation error, not as medium/GEM infeasibility.

## Step 6: final result

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

head(result$reaction_ranking)
head(result$condition_contrast)
```

Final files: `step_comparison.rds` and `regcompass_result.rds`.

## Error interpretation

- **`RNA normalized assay data are not aligned`**: column-order differences are accepted in 1.8.2. A remaining error means the cell-ID sets genuinely differ.
- **`Some features contain 0 total counts`**: cell-type-local all-zero peaks are now omitted from TF-IDF and restored as zeros.
- **`The medium-constrained parent GEM is not feasible: error`**: reinstall 1.8.2 and confirm `highs` is installed. Solver availability is checked before feasibility analysis.
