# RegCompassR

RegCompassR integrates paired single-cell RNA+ATAC regulatory evidence with a
shared metabolic model for reaction-level analysis.

## Automatic analysis modes

Stage 1 selects the GRN mode from the supplied condition metadata:

| Input condition design | GRN implementation | Condition coefficients |
|---|---|---|
| Two or more levels | `Pando::infer_condition_grn()` | Yes |
| One level | original `Pando::infer_grn()` | No |
| Column omitted or absent | original `Pando::infer_grn()` | No |

The standard-Pando route is not an RNA-only fallback. Standard TF–peak–gene
coefficients are projected as cell-level TF RNA × peak ATAC terms, aggregated by
exact SuperCell membership, and then passed through the same RNA-support, GPR,
reaction-expression, penalty, and metabolic-LP stages. No artificial condition
coefficient or contrast is created.

With two or more conditions, the canonical unversioned
`pando_condition_grn_fit` contract supplies outer-heldout common-support
projections. Conditions are represented by absolute coefficients on one
within-cell-type equal-condition coordinate.

## Cell-type-independent, condition-joint SuperCell construction

RegCompass separates graph scope from metacell purity. It does not create a
temporary `condition__cell_type` metadata field and does not build a graph per
condition. Instead, RNA and ATAC embedding blocks are standardized within each
broad cell type using all conditions of that cell type, then passed to:

```r
SuperCell::SCimplify_by_graph_group_from_embedding(
  X = joint_rna_pca_atac_lsi_embedding,
  cell.graph.group = cell_type,
  cell.split.condition = condition,
  gamma = 30
)
```

The invariants are:

- each broad cell type has an independent kNN graph;
- all conditions of that cell type are jointly present during embedding
  standardization, neighbour search, and graph clustering;
- condition is applied after clustering to split mixed preliminary memberships;
- final metacells are pure for both cell type and condition;
- no `sample` column is used for graph construction or membership splitting.

For an omitted condition, `cell.split.condition` is `NULL`. RNA and ATAC raw
counts are then aggregated from the returned `cell_id → metacell_id` membership.
The cache and stage contract records:

```text
graph_scope = one_independent_graph_per_cell_type
condition_scope = all_conditions_joint_within_cell_type_graph
membership_split_timing = after_joint_graph_clustering
embedding_scaling = within_celltype_joint_condition_equal_modality_blocks
temporary_combined_stratum = FALSE
```

## Workflow

```text
paired RNA+ATAC cells
→ automatic standard/condition-aware Pando routing
→ one multimodal graph per cell type with conditions jointly embedded
→ condition-pure metacell membership after graph clustering
→ cell-type Gamma–Poisson RNA latent expression
→ cell-first TF×ATAC regulatory projection
→ GPR reaction expression and COMPASS-like penalties
→ shared medium-specific GEM and directional LP scoring
```

Canonical behavior is implemented directly in the defining functions. The
package no longer depends on `zzz` runtime function replacement or schema-body
rewriting.

## Installation

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@agent/celltype-joint-condition-graphs")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass@agent/celltype-joint-condition-graphs")
```

The branch-qualified SuperCell dependency is temporary while the companion PR is
under review. After that PR is merged, install the merged `Supercell2` branch.

## Required input

A paired-cell Seurat object with:

- RNA and ATAC count assays for the same cells;
- a broad-cell-type metadata column;
- optional condition metadata;
- RNA PCA and ATAC LSI reductions;
- genome-compatible peak coordinates.

## One-shot run with multiple conditions

```r
library(RegCompassR)
library(BSgenome.Hsapiens.UCSC.hg38)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

result <- rc_run_regcompass(
  object = A,
  gem = gem,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 100L,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      condition_weight = "equal",
      outer_nfolds = 5L,
      inner_nfolds = 5L,
      scale = TRUE
    )
  ),
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L
  ),
  layer1_args = list(
    comparison_support = "auto",
    regulatory_alpha = 1,
    gpr_and_method = "min"
  )
)
```

## One-shot run without a multi-level condition

Either omit the condition column:

```r
result <- rc_run_regcompass(
  object = A,
  gem = gem,
  outdir = "RegCompass_standard",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = NULL,
  celltype_col = "cell_type"
)
```

or supply a column containing one level. Both cases automatically use original
Pando `infer_grn()`. Standard-only Pando arguments can be supplied through
`pando_args$pando_infer_args`; condition-only arguments are not required.

## Main outputs

```r
result$analysis_mode
result$condition_coefficients_calculated
result$grn
result$metacells$membership
result$metacells$input_design
result$layer1$projection_provenance
result$microcompass$model_cache_summary
result$reaction_ranking
result$condition_summary
result$condition_contrast
```

For one effective condition, `condition_contrast` is empty and
`reaction_ranking` remains available. Metacells are within-dataset statistical
units, not biological replicates.

## Documentation

- [Workflow](docs/workflow.md)
- [Stage contracts](docs/stage-interface-contracts.md)
- [Metacell graph contract](docs/metacell-graph-contract.md)
- [Run modes](docs/run-modes-and-stepwise-workflow.md)
- [Mathematical model](docs/mathematical-model.md)
- [Public API](docs/functions.md)
