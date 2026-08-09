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
│  └─► Infer cell-type regulatory networks; use common-dictionary condition GRNs
│      when at least two conditions are retained, otherwise use standard Pando.
│
├─ 2. Multimodal metacells
│  └─► Aggregate paired cells within each cell type into multimodal metacells while
│      preserving condition purity for downstream comparison.
│
├─ 3. GPR-supported reaction meta-modules
│  └─► Convert active GRN-supported metabolic genes into complete-GPR core reactions
│      and expand them by subsystem and direct reaction-database equivalence.
│
├─ 4. Reaction evidence projection
│  └─► Project RNA-only and RNA+ATAC regulatory evidence onto GEM reactions to obtain
│      directly comparable reaction-level evidence for each metacell.
│
├─ 5. Structural model + directional scoring
│  └─► Build cell-type-specific structural metabolic models and calculate directional
│      penalties while reusing the same model for the matched RNA-only control.
│
└─ 6. Result assembly and condition comparison
   └─► Assemble reaction evidence, rankings and condition contrasts for biological
       interpretation and downstream analysis.
```

Stage-specific behavior:

- **GRNs:** common-dictionary condition fits require at least two retained conditions within a cell type. RegCompass parallelizes pooled-background and condition × cell-type candidate discovery, waits for all candidate tasks of each cell type, freezes one exact edge dictionary per cell type, and then parallelizes condition × cell-type fixed-dictionary GLMs. Standard Pando is parallelized across broad cell types and uses the fixed `tf_cor` and `peak_cor` values supplied through `pando_infer_args` without any cell-count-dependent threshold adjustment.
- **Meta-modules:** complete-GPR core reactions are expanded by core subsystem and direct KEGG, Reactome or master-Rhea equivalence.
- **Structural models:** `meta_module_gem` builds one model per cell-type × medium combination; CORDA2 is the default completion route and runs without a structural time limit. `model_params$completion_time_limit` is rejected for CORDA2 and remains available only for supplementary non-CORDA2 completion such as FASTCORE. `full_gem` keeps the complete GEM with medium bounds.
- **Matched control:** RNA+ATAC and RNA-only penalties reuse the same structural model, bounds, medium and target directions.

Stage 3 retains only newly derived meta-module information and does not serialize the Stage 1 GRN payload again.

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
      tf_cor = 0.1,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  workers = 10L
)
```

`workers` is the only workflow-level parallel setting. The default is `10L`. RegCompass automatically selects `SnowParam(type = "SOCK")` on Windows and `MulticoreParam` on Linux/macOS; the effective cap is `min(workers, max(1, detected logical CPUs - 2))` and each dispatch shrinks further to its number of independent jobs.

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
