# RegCompassR

RegCompassR 1.8.1 runs the following RNA+ATAC workflow:

```text
single-cell RNA normalization
→ ATAC TF-IDF shared across conditions within each cell type
→ Pando GRN for each condition × cell type
→ condition-only SuperCell2 metacells
→ GRN-derived core reactions and meta-modules
→ RNA+ATAC reaction expression
→ directional COMPASS-like scoring
```

The canonical defaults are `peak_cor = 0.01` for Pando and `gamma = 75` for SuperCell2. Sample metadata are optional and are not used for balancing, weighting, downsampling, or grouping.

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

A locally downloaded Pando source package is also supported:

```r
install.packages(
  "~/Pando_regcompass.tar.gz",
  repos = NULL,
  type = "source"
)
```

RegCompass validates the required Pando API. GitHub remote metadata are not required for a local or offline source installation.

## Required input

`object` must be a paired-cell Seurat multiome object with:

- an RNA assay containing raw counts;
- an ATAC `ChromatinAssay` containing peak counts for the same cell IDs;
- complete condition and cell-type metadata;
- peak coordinates and genome build matching `genome`;
- a PFM/PWM collection accepted by Pando/motifmatchr. Use `Pando::motifs`; do not pass the `motif2tf` annotation table as `pfm`.

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

stopifnot(
  inherits(A, "Seurat"),
  all(c("RNA", "ATAC") %in% names(A@assays)),
  inherits(A[["ATAC"]], "ChromatinAssay"),
  all(c("dataset", "epithelial_or_stem") %in% colnames(A@meta.data)),
  !anyNA(A$dataset),
  !anyNA(A$epithelial_or_stem)
)

gem <- rc_prepare_gem(species = "human", version = "2.0.0")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "high_glucose",
  species = "human"
)
```

RNA and ATAC normalized matrices are aligned by cell name; different column order is accepted. Peaks absent from one cell type are retained as exact zeros but are excluded from that cell type's TF-IDF calculation.

## One-shot run

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
  gem = gem,
  medium_scenarios = medium_scenarios,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  pando_args = list(
    min_cells = 100,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr"
    )
  ),
  metacell_args = list(
    gamma = 75,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  layer1_args = list(local_fastcore = TRUE),
  layer2_args = list(target_direction = "both", solver = "highs")
)
```

The default `highs` backend is a required dependency. If another solver is selected, RegCompass checks its R package before constructing the medium-constrained model and reports a solver-installation error separately from biological infeasibility.

## Main outputs

- `01_single_cell_grn/pando_group_status.tsv.gz`: one row per condition × cell-type GRN.
- `01_single_cell_grn/pando_tf_peak_gene_significant.tsv.gz`: significant Pando edges.
- `02_condition_metacells/metacell_metadata.tsv.gz`: final metacell labels and purity diagnostics.
- `03_meta_modules/grn_metacell_group_coverage.tsv.gz`: GRN-to-metacell coverage validation.
- `03_meta_modules/core_gene_reaction.tsv.gz`: complete-GPR core reactions.
- `03_meta_modules/meta_module_reactions.tsv.gz`: biological meta-module membership before FASTCORE support.
- `04_layer1/step_layer1.rds`: RNA support, ATAC modifier, integrated gene support, and reaction expression.
- `05_layer2/step_layer2.rds`: directional LP scores, penalties, feasibility, and diagnostics.
- `06_results/regcompass_result.rds`: assembled final result.

See [the stepwise tutorial](docs/run-modes-and-stepwise-workflow.md) for stage-by-stage execution and inspection.
