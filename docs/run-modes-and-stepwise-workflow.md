# RegCompassR stepwise workflow

This tutorial describes the canonical 1.8.2 workflow and the files that connect each stage.

## 1. Installation

```r
install.packages(c("remotes", "highs"))
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0", upgrade = "never")
remotes::install_github("1667857557/Pando_regcompass", upgrade = "never")
remotes::install_github("1667857557/Regcompass", upgrade = "never")
```

Pando may instead be installed from a local source archive:

```r
install.packages("~/Pando_regcompass.tar.gz", repos = NULL, type = "source")
```

A local installation does not contain GitHub remote metadata. RegCompass checks the required Pando functions and accepts the package without a remote-origin warning.

## 2. Input contract

The same single-cell Seurat object is supplied to Steps 1 and 2. It must contain paired RNA and ATAC measurements for the same cell identities.

Required components:

| Input | Requirement |
|---|---|
| RNA assay | Raw counts with gene symbols compatible with the selected GEM |
| ATAC assay | `ChromatinAssay` peak-count matrix |
| Cells | RNA, ATAC, and metadata refer to the same cell IDs; order may differ |
| Condition | Complete metadata column, for example `dataset` |
| Cell type | Complete metadata column, for example `epithelial_or_stem` |
| Genome | Matches ATAC coordinates, for example hg38 |
| `pfm` | PFM/PWM collection accepted by Pando/motifmatchr; `Pando::motifs`, not `motif2tf` |
| Fragments | Not required when `fragment_files = FALSE`; existing peak counts are aggregated |

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
```

RegCompass runs `NormalizeData` on RNA. For ATAC, it calculates one TF-IDF reference per cell type using all conditions. A peak that is globally present but absent in one cell type remains zero for that cell type and is not passed to `RunTFIDF`, preventing zero-total warnings.

## 3. GEM and medium

```r
gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "high_glucose",
  species = "human"
)

rc_validate_gem(gem)
head(medium_scenarios)
```

The default solver is `highs`, which is a package dependency. RegCompass checks the selected solver before constructing a medium-constrained GEM. A missing solver is therefore reported as an installation error, not as model infeasibility.

## 4. Step 1: single-cell GRNs

**Input:** original Seurat object, GEM, motifs, genome, condition and cell-type columns.

**Computation:** global RNA normalization, cell-type-shared ATAC TF-IDF, then one Pando model per `condition × cell type`.

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

Primary files:

- `single_cell_grn.rds`: GRN result used by Step 3.
- `step_grn.rds`: restartable stage object.
- `pando_group_status.tsv.gz`: group cell counts and edge counts.
- `pando_tf_peak_gene_all.tsv.gz`: all fitted coefficients.
- `pando_tf_peak_gene_significant.tsv.gz`: filtered edges used downstream.
- `pando_objects/*.rds`: optional fitted Pando objects.

## 5. Step 2: condition-only metacells

**Input:** the original, unmodified Seurat object.

**Computation:** SuperCell2 is stratified only by condition. Cell type is assigned afterwards from membership. An exact tie between two cell types is rejected because no unique GRN can be selected.

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

with(
  step2$pooled$metacell_meta,
  table(.data[[condition_col]], .data[[celltype_col]])
)
summary(step2$pooled$metacell_meta$dominant_celltype_fraction)
```

For base R compatibility, the table can also be written as:

```r
with(
  step2$pooled$metacell_meta,
  table(
    step2$pooled$metacell_meta[[condition_col]],
    step2$pooled$metacell_meta[[celltype_col]]
  )
)
```

Primary files:

- `step_metacells.rds`: restartable stage object.
- `merged_metacell_object.rds`: normalized RNA+ATAC metacell Seurat object.
- `metacell_metadata.tsv.gz`: final condition, dominant cell type, purity, and mixed-cell diagnostics.
- `metacell_membership.tsv.gz`: single-cell-to-metacell assignments.
- `metacell_celltype_composition.tsv.gz`: full cell-type composition per metacell.

## 6. Step 3: core reactions and meta-modules

**Input:** Step 1 GRNs, Step 2 metacells, and the same GEM.

**Computation:** validates bidirectional GRN/metacell group coverage, maps GRN metabolic genes to complete-GPR core reactions, expands through subsystem and shared KEGG/Reactome/master-Rhea identifiers, then adds local FASTCORE support reactions.

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

Primary files:

- `grn_metacell_group_coverage.tsv.gz`: required GRN/metacell alignment.
- `core_gene_reaction.tsv.gz`: complete-GPR core mappings.
- `meta_module_reactions.tsv.gz`: biological membership before FASTCORE support.
- `condition_meta_modules.rds`: condition × cell-type modules.
- `global_meta_modules.rds`: deduplicated union used by Layer 2.
- `local_fastcore/`: completion models and diagnostics.

## 7. Step 4: RNA+ATAC reaction expression

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  tau = 0.20
)

dim(step4$gene_support_rna)
dim(step4$gene_regulatory_modifier)
dim(step4$reaction_expression)
stopifnot(identical(
  colnames(step4$reaction_expression),
  step4$unit_meta$pool_id
))
```

`step_layer1.rds` contains RNA support, ATAC regulatory modifiers, integrated gene support, parsed GPRs, reaction expression, and metacell metadata.

## 8. Step 5: directional scoring

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

Primary outputs include `score`, `penalty`, `vmax`, `feasible`, LP diagnostics, model diagnostics, and the persistent model cache.

## 9. Step 6: final result

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

names(result)
head(result$reaction_ranking)
head(result$condition_contrast)
```

The final files are `step_comparison.rds` and `regcompass_result.rds`.

## Common errors

### `RNA normalized assay data are not aligned`

Different column order is valid and is now automatically corrected. A remaining error means RNA and the Seurat object contain genuinely different cell IDs.

### `Some features contain 0 total counts`

Globally nonzero peaks can be absent from one cell type. RegCompass now excludes those local all-zero rows from that TF-IDF calculation and restores them as exact zeros.

### `The medium-constrained parent GEM is not feasible: error`

An `error` status may indicate a missing solver rather than biological infeasibility. Version 1.8.2 requires `highs` and performs a solver preflight. True infeasibility is reported only after a solver has run successfully.
