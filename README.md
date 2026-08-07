# RegCompassR

RegCompassR combines paired single-cell RNA and ATAC regulatory evidence with genome-scale metabolic models to calculate cell-type-resolved reaction scores. The outputs are model-derived priorities and penalties, not direct measurements of intracellular flux.

## Installation

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

## Required input

Use a paired-cell Seurat object containing:

- RNA and ATAC count assays for the same cells;
- a broad cell-type metadata column;
- RNA PCA or Harmony and ATAC LSI reductions;
- genome-compatible peak coordinates;
- an optional condition column.

Stage 1 retains analysis groups meeting `pando_args$min_cells` (300 by default). Within each broad cell type, two or more retained conditions use common-dictionary condition GRNs; one retained condition or no condition uses standard Pando.

## Workflow

```text
Paired RNA + ATAC cells, GEM and medium
                  │
                  ▼
1. Cell-type GRN inference
                  │
                  ▼
2. Condition-pure multimodal metacells
                  │
                  ▼
3. GPR-supported reaction meta-modules
                  │
                  ▼
4. RNA-only and RNA+ATAC reaction evidence
                  │
                  ▼
5. Structural model construction and directional scoring
                  │
                  ▼
6. Reaction tables, rankings and condition comparisons
```

Stage-specific behavior:

- **GRNs:** common-dictionary condition fits are used only when a cell type contains at least two retained conditions; otherwise the workflow uses standard Pando.
- **Meta-modules:** core reactions require complete GPR support from active metabolic target genes. Expansion uses core subsystems and direct KEGG, Reactome or master-Rhea reaction equivalence.
- **Structural models:** `meta_module_gem` constructs a model for each cell-type × medium combination. CORDA2 is the default completion route and FASTCORE is optional. `full_gem` retains the complete GEM and applies medium bounds without context-specific reconstruction.
- **Matched control:** the primary RNA+ATAC penalty and RNA-only control reuse the same structural model, bounds, medium and target directions. Their difference describes the effect of including regulatory evidence within this model; it should not be interpreted as causal isolation or experimental flux validation.

Runtime RDS artifacts retain the `.rds` extension and are written with gzip compression. Stage 1 stores the complete checkpoint as `step_grn.rds`; the redundant `single_cell_grn.rds` intermediate is not written.

## Simple workflow

```r
library(RegCompassR)

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300L,
    pando_infer_args = list(
      tf_cor = 0.05,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  )
)
```

When condition metadata are unavailable, omit `condition_col` or set `condition_col = NULL`. The workflow then uses standard Pando with one internal background label.

The one-shot workflow prepares the species GEM and default plasma-like medium, runs the six stages above and returns condition-level comparisons only when the design contains multiple retained conditions.

Main outputs:

```r
result$grn$cell_type_analysis_mode
result$reaction_catalog
result$reaction_evidence
result$reaction_ranking
result$condition_contrast
result$reaction_comparison_by_metacell
```

## Documentation

1. [Quick run](docs/tutorial-01-quick-start.md)
2. [Stepwise run](docs/tutorial-02-stepwise-audit.md)
3. [Post analysis](docs/tutorial-04-post-analysis.md)
4. [Principles and mathematical formulas](docs/mathematical-model.md)
5. [Function reference](docs/functions.md)
