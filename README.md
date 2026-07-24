# RegCompassR

RegCompassR 1.8.2 implements a GRN-first RNA+ATAC metabolic workflow:

```text
single-cell RNA normalization
→ cell-type-shared ATAC TF-IDF
→ Pando GRNs for each condition × cell type
→ condition-level, label-guided SuperCell2 metacells
→ complete-GPR core reactions and annotation-expanded meta-modules
→ local FASTCORE completion and one global union GEM
→ RNA+ATAC reaction expression
→ directional COMPASS-like LP scoring
→ optional direct database-linked non-core scoring in the same union GEM
```

Canonical defaults are `peak_cor = 0.01` and `gamma = 75`. Sample metadata are provenance only; they are not used for balancing, weighting, or grouping.

## Installation

```r
install.packages(c("remotes", "highs", "BiocManager", "ggplot2"))
BiocManager::install("BiocParallel", ask = FALSE, update = FALSE)
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0", upgrade = "never")
remotes::install_github("1667857557/Pando_regcompass", upgrade = "never")
remotes::install_github("1667857557/Regcompass", upgrade = "never")
```

## Input contract

`object` must contain paired RNA and ATAC counts for the same cells, complete condition and cell-type metadata, a Signac `ChromatinAssay`, and peak coordinates matching `genome`. Use `Pando::motifs` as `pfm`; `motif2tf` is not a motif matrix collection.

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

gem <- rc_prepare_gem(species = "human", version = "2.0.0")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "physiologic",
  species = "human"
)
```

`physiologic` resolves to the species-specific plasma preset. Culture, nutrient-sensitivity, technical, and custom media are documented in [Predefined extracellular medium scenarios](docs/medium-presets.md).

## One-shot Linux run

Set numerical-library threads to one before starting R when using multiple workers.

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

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
      adjust_method = "fdr",
      parallel = FALSE
    )
  ),
  metacell_args = list(
    gamma = 75,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  layer1_args = list(
    local_fastcore = TRUE,
    local_fastcore_args = list(solver = "highs", parallel = TRUE)
  ),
  layer2_args = list(target_direction = "both", solver = "highs"),
  upstream_workers = 16L,
  layer2_workers = 12L,
  parallel_backend = "multicore"
)
```

The workflow validates the solver, stage classes, GEM fingerprint, workflow metadata, and ordered metacell IDs before connecting stages.

## Score direct database links of selected cores

After a stepwise `meta_module_gem` Layer 2 run, select previous core reactions or their GPR genes. The selected cores are mapping anchors only. The function directly identifies reactions sharing a KEGG, Reactome, or master-Rhea ID with each selected core and runs a second LP only for linked reactions that were not global core targets in the original Layer 2 run.

Same-subsystem expansion, recursive propagation through intermediate reactions, FASTCORE-only support, generic union members, metabolite-neighbour reactions, and previously scored global cores are excluded.

```r
expanded <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_expanded_target_scoring",
  core_genes = c("GCLC", "GCLM", "GSS", "GSR", "G6PD", "PGD"),
  gene_match = "complete_gpr"
)
```

`expanded$expanded_reaction_catalog` records `anchor_core_reaction_id`, mapping type, database identifier, and core-exclusion status. `expanded$expanded_scoring_targets` contains unique direct database-linked non-core targets. `expanded$microcompass$penalty` is the primary result; lower penalty means stronger evidence-supported flux compatibility.

## Compare the same reaction across conditions

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  cell_types = "stem-cell_like",
  p_adjust_scope = "celltype_contrast_medium"
)

condition_stats$omnibus
condition_stats$pairwise

p <- rc_plot_condition_reaction(
  result,
  reaction_id = "MAR06231",
  cell_type = "stem-cell_like",
  target_direction = "reverse",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  annotation_p = "p_adj"
)
```

The plot shows one point per metacell and adjusted significance brackets. These are metacell-level, within-dataset comparisons rather than biological-replicate inference. See [Condition-associated reaction statistics](docs/condition-reaction-statistics.md).

## Main outputs

- `01_single_cell_grn/`: GRN status and Pando edges.
- `02_condition_metacells/`: metacell counts, membership, labels, and purity.
- `03_meta_modules/`: core reactions, annotation-expanded modules, and FASTCORE diagnostics.
- `04_layer1/step_layer1.rds`: RNA support, ATAC modifier, GPR diagnostics, and reaction expression.
- `05_layer2/step_layer2.rds`: penalties, relative scores, `vmax`, feasibility, model cache, and LP diagnostics.
- `06_results/regcompass_result.rds`: rankings, annotations, and evidence provenance.
- optional target-union output: selected core anchors, direct database relation catalog, unique non-core scoring targets, source-model hashes, and second-pass LP results.

## Tutorials

| Level | Use | Tutorial |
|---|---|---|
| 1 | minimal validated one-shot run | [Quick start](docs/tutorial-01-quick-start.md) |
| 2 | stage-by-stage run and audit gates | [Stepwise audit](docs/tutorial-02-stepwise-audit.md) |
| 3 | restart, sensitivity, resources, and failure diagnosis | [Advanced restart](docs/tutorial-03-advanced-restart.md) |
