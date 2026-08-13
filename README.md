# RegCompassR

RegCompassR integrates paired single-cell RNA and ATAC regulatory evidence with genome-scale metabolic models and returns cell-type-resolved reaction scores. The scores are model-derived reaction support/penalty measures, not direct flux measurements.

## Installation

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

## Required input

Use a paired-cell Seurat object with RNA and ATAC count assays for the same cells, broad cell-type metadata, RNA PCA or Harmony, ATAC LSI, and genome-compatible peak coordinates. A condition column is optional.

## Workflow

```text
Paired RNA + ATAC cells + GEM + medium
│
├─ 1. GRN inference
│  └─► Infer cell-type regulatory networks. Cell types with at least two retained
│      conditions use one correlation-screened common dictionary with condition-
│      specific no-fusion ridge fits; one-condition cell types use standard Pando.
│
├─ 2. Multimodal metacells
│  └─► Aggregate paired cells within each broad cell type into multimodal metacells
│      while preserving condition-pure membership for downstream comparison.
│
├─ 3. GPR-supported reaction meta-modules
│  └─► Map regulatory-supported metabolic genes through complete GPR rules to core
│      reactions and expand the reaction context with the package meta-module rules.
│
├─ 4. Reaction evidence projection
│  └─► Build quantitative RNA+regulatory reaction evidence and the matched RNA-only
│      control on the same metacell/reaction coordinate system.
│
├─ 5. Structural model + directional scoring
│  └─► Build cell-type × medium structural GEMs with CORDA2, preserve required core
│      reactions, compute directional Vmax on the final GEM, and solve Step 2 penalties.
│      The matched RNA-only control reuses the same structural model and cached Vmax.
│
└─ 6. Result assembly and condition comparison
   └─► Assemble reaction evidence, directional scores, rankings, annotations and
       condition contrasts for downstream biological interpretation.
```

Current workflow contracts:

- **GRNs:** all conditions of a multi-condition cell type are fitted on the same exact common edge dictionary, so condition-specific coefficients are directly aligned by edge. Edge activity is determined by that condition's own estimable BH-supported ridge evidence; global/local correlation support remains dictionary provenance.
- **Reaction evidence:** the quantitative COMPASS-like penalty path and the bounded structural-confidence path are kept distinct. RNA-only is retained as a matched control rather than a separate structural model.
- **Structural scoring:** CORDA2 reconstructs the cell-type × medium model first. Directional Vmax is evaluated on that final GEM and reused by Step 2; RNA-only scoring changes the evidence-dependent objective rather than rebuilding the model or recomputing Vmax.
- **Parallelism:** one top-level `workers` cap is shared across workflow stages; heavy reaction-level work uses the available worker budget without nested solver threading.

## Minimal workflow

```r
library(RegCompassR)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  medium_scenarios = medium_scenarios,
  condition_col = "condition",
  celltype_col = "cell_type",
  workers = 10L
)
```

For condition-aware analysis, cell types with at least two retained conditions use the condition-GRN route; cell types with one retained condition use standard Pando. `condition_col = NULL` is valid when condition metadata are absent.

## Main defaults

- Human GEM: Human-GEM `2.0.0`; mouse GEM: Mouse-GEM `1.8.0`.
- Human default medium: `normal_human_plasma`; mouse default medium: `mouse_plasma`.
- Stage 1 retained-group threshold: `pando_args$min_cells = 500L`.
- Stage 2 public retained-stratum threshold: `metacell_args$min_cells_per_stratum = 500L`.
- Structural mode: `model_mode = "meta_module_gem"`, with CORDA2 as the default completion route.
- Parallelism: one top-level `workers` cap, default `10L`.

These thresholds are configurable through the current public arguments.

## Medium presets

`rc_make_medium_scenarios()` supports `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`.

See [`docs/medium-presets.md`](docs/medium-presets.md) for compositions, provenance, and custom-medium formats.

## Documentation

- [Quick start](docs/tutorial-01-quick-start.md)
- [Restartable stepwise workflow](docs/tutorial-02-stepwise-audit.md)
- [Post analysis](docs/tutorial-04-post-analysis.md)
- [Function reference](docs/functions.md)
- [Mathematical specification](docs/mathematical-model.md)

Tutorials and Rd pages document interfaces. Equations and quantitative definitions are maintained only in `docs/mathematical-model.md`.
