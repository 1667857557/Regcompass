# RegCompassR

RegCompassR combines paired single-cell RNA and ATAC regulatory evidence with genome-scale metabolic models to calculate cell-type-resolved reaction scores. Outputs are model-derived priorities and penalties, not direct intracellular flux measurements.

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

Stage 1 keeps groups meeting `pando_args$min_cells` (300 by default). A cell type with at least two retained conditions uses common-dictionary condition GRNs; otherwise it uses standard Pando.

## Workflow

```text
Paired RNA + ATAC cells + GEM + medium
│
├─ 1. GRN inference
│  ├─ ≥2 retained conditions within a cell type ──► common-dictionary condition GRNs
│  └─ 0–1 retained condition within a cell type ──► standard Pando GRN
│                                                   └─► step_grn.rds
│
├─ 2. Multimodal metacells
│  └─ one multimodal graph per cell type; memberships remain condition-pure
│                                                   └─► step_metacells.rds
│                                                       merged_metacell_object.rds
│
├─ 3. GPR-supported reaction meta-modules
│  ├─ active metabolic targets ─► complete-GPR core reactions
│  └─ core reactions ───────────► subsystem + KEGG/Reactome + master-Rhea expansion
│                                                   ├─► condition_meta_modules.rds
│                                                   │    Stage 3 information only;
│                                                   │    Stage 1 GRN is not copied
│                                                   ├─► merged_meta_modules.rds
│                                                   └─► step_meta_modules.rds
│
├─ 4. Reaction evidence projection
│  ├─ RNA expression ───────────► RNA-only reaction evidence
│  └─ RNA + Pando regulation ───► RNA+ATAC reaction evidence
│                                                   └─► step_layer1.rds
│
├─ 5. Structural model + directional scoring
│  ├─ meta_module_gem ──────────► one structural model per cell type × medium
│  ├─ full_gem ─────────────────► full GEM with medium bounds
│  └─ matched RNA-only control ─► same model, bounds and target directions
│                                                   └─► step_layer2.rds
│                                                       model_cache/
│
└─ 6. Result assembly and condition comparison
   ├─ reaction catalogue / evidence / rankings
   ├─ condition summaries and contrasts
   └─ metacell-level RNA+ATAC vs RNA-only comparison
                                                    ├─► regcompass_result.rds
                                                    ├─► step_comparison.rds
                                                    └─► reaction_*.tsv.gz
```

Stage-specific behavior:

- **GRNs:** common-dictionary condition fits require at least two retained conditions within a cell type; other cell types use standard Pando.
- **Meta-modules:** complete-GPR core reactions are expanded by core subsystem and direct KEGG, Reactome or master-Rhea equivalence.
- **Structural models:** `meta_module_gem` builds one model per cell-type × medium combination; CORDA2 is the default completion route and FASTCORE is optional. `full_gem` keeps the complete GEM with medium bounds.
- **Matched control:** RNA+ATAC and RNA-only penalties reuse the same structural model, bounds, medium and target directions.

Runtime RDS artifacts use gzip compression. Stage 1 stores the GRN checkpoint in `step_grn.rds`. Stage 3 stores only newly derived meta-module information in `condition_meta_modules.rds`; it does not serialize the Stage 1 GRN again.

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

When condition metadata are unavailable, omit `condition_col` or set `condition_col = NULL`; the workflow then uses standard Pando with one internal background label.

The one-shot workflow prepares the species GEM and default plasma-like medium, runs the six stages above and returns condition comparisons only when multiple retained conditions are available.

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
