# RegCompassR

RegCompassR 1.8.1 uses a GRN-first RNA+ATAC workflow:

```text
single-cell RNA NormalizeData across all cells
→ ATAC TF-IDF within each cell type across conditions
→ one Pando GRN per condition × cell type (peak_cor = 0.01 by default)
→ condition-only SuperCell2 metacells (gamma = 75 by default)
→ unambiguous post hoc dominant cell-type labels from metacell membership
→ validate GRN ↔ metacell condition × cell-type coverage
→ complete-GPR core reactions from each condition × cell-type GRN
→ subsystem + KEGG/Reactome + master-Rhea expansion
→ local FASTCORE feasibility completion
→ RNA+ATAC reaction expression
→ directional COMPASS-like minimum-penalty scoring
```

Sample metadata are optional. They are not used for sample balancing, downsampling, weighting, or metacell grouping. Cell type is also not used to stratify metacell construction; it is assigned afterwards from the dominant member-cell label. Purity and mixed-cell-type diagnostics are retained, and exact dominant-cell-type ties are rejected because no condition × cell-type GRN can be assigned unambiguously.

## Installation

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0", upgrade = "never")
remotes::install_github("1667857557/Pando_regcompass", upgrade = "never")
remotes::install_github("1667857557/Regcompass", upgrade = "never")
```

## One-shot workflow

```r
result <- rc_run_regcompass_one_shot(
  object = object,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
  condition_col = "condition",
  celltype_col = "cell_type",
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
  layer1_args = list(
    regulatory_alpha = 1,
    tau = 0.20,
    local_fastcore = TRUE
  ),
  layer2_args = list(target_direction = "both", solver = "highs")
)
```

## Stepwise workflow

```r
step1 <- rc_regcompass_step_grn(
  object = object,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "condition",
  celltype_col = "cell_type",
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

step2 <- rc_regcompass_step_metacells(
  object = object,
  outdir = "RegCompass_steps/02_metacells",
  condition_col = "condition",
  celltype_col = "cell_type",
  fragment_files = FALSE,
  metacell_args = list(gamma = 75)
)

step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules",
  layer1_args = list(local_fastcore = TRUE)
)

step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1"
)

step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_layer2"
)

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

## Auditable stage outputs

The restartable stages write the following primary contracts:

- `01_grn`: `single_cell_grn.rds`, `step_grn.rds`, Pando group status and all/significant edge tables, plus optional per-group Pando objects.
- `02_metacells`: `step_metacells.rds`, `merged_metacell_object.rds`, final `metacell_metadata.tsv.gz`, `metacell_membership.tsv.gz`, and cell-type composition/summary tables. These root-level tables contain the post hoc labels used downstream, unlike the raw per-stratum intermediates.
- `03_meta_modules`: `grn_metacell_group_coverage.tsv.gz`, condition-specific and global meta-module RDS files, core-reaction and reaction-membership tables, and local FASTCORE diagnostics.
- `04_layer1`: `step_layer1.rds`, containing metacell RNA support, ATAC regulatory modifiers, integrated gene support, and reaction expression.
- `05_layer2`: exported microCOMPASS matrices/tables, `step_layer2.rds`, and a persistent model cache.
- `06_results`: `step_comparison.rds` and `regcompass_result.rds`, retaining both condition-specific and global meta-modules.

The Pando coefficients are learned from single cells. Metacells are built only within condition and receive dominant cell-type labels afterwards so the corresponding condition × cell-type GRN can be applied in Layer 1. Before meta-module construction, every successful GRN group must have at least one scoring metacell and every scoring metacell group must have a successful GRN with significant edges. Meta-module expansion remains restricted to complete-GPR core reactions, their subsystems, shared KEGG or Reactome identifiers, and identical master-Rhea identifiers. Local FASTCORE adds only feasibility-support reactions.
