# RegCompassR

RegCompassR provides one supported RNA+ATAC workflow:

```text
condition × sample × cell type
→ metacells and stratum-specific Pando GRNs
→ local FASTCORE meta-modules
→ sample-balanced global calibration and one shared GEM
→ directional microCOMPASS scoring
```

## Installation

```r
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

## Quick start

The one-shot entry point prepares Human-GEM and the shared model-bound medium
when they are not supplied:

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)
data(motifs, package = "Pando")

result <- rc_run_regcompass_one_shot(
  object = object,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  metacell_args = list(
    gamma = 50,
    min_cells_per_stratum = 1000,
    min_metacell_size = 10
  ),
  pando_args = list(
    min_metacells = 20,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr"
    )
  ),
  layer1_args = list(
    local_fastcore = TRUE,
    sample_balance = TRUE,
    expression_batch_correction = "none"
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  )
)
```

## Choosing analysis parameters

The values above are starting points, not fixed analysis parameters. Choose
them before a full run from stratum sizes, expected biological heterogeneity,
and a small pilot run.

| Setting | Role | Selection guidance |
| --- | --- | --- |
| `metacell_args$gamma` | Approximate cells represented by each metacell | Use a smaller value for more metacells and finer heterogeneity, or a larger value for stronger aggregation. Keep one value across all strata so resolution is comparable. |
| `metacell_args$min_cells_per_stratum` | Minimum cells in each condition × sample × cell-type stratum | Set high enough to avoid unstable strata, but check that every biological sample retains at least one stratum. |
| `metacell_args$min_metacell_size` | Flags undersized, low-power metacells | Increase when sparse RNA/ATAC profiles are unstable; do not use it to compensate for an unsuitable `gamma`. |
| `pando_args$min_metacells` | Minimum metacells required for Pando | It must be compatible with `floor(n_cells / gamma)`. Strata below it are skipped, so inspect `00_strata/stratum_workflow_status.tsv.gz` after a pilot. |
| `pando_infer_args` | Pando model and correlation/FDR filters | Start with the shown GLM settings; tighten correlation thresholds only when the pilot produces excessive weak edges. Apply one policy to every stratum. |
| `layer1_args$local_fastcore` | Completes each local metabolic module | Keep enabled for the canonical path. `sample_balance = TRUE` prevents samples with more metacells from dominating global calibration. |
| `layer1_args$expression_batch_correction` | Optional technical-batch correction | Keep `"none"` unless a documented technical batch exists. If using `"limma"`, provide technical and preserved biological design columns; never remove `sample_id` as batch. |
| `layer2_args$target_direction` | Forward, reverse, or both-direction scoring | Use `"both"` unless the GEM direction or scientific question justifies one direction. |
| `layer2_args$solver` | LP solver | `"highs"` is the default open-source choice; use Gurobi only in a licensed environment. |

A practical pilot is to tabulate cells per strict stratum, choose one `gamma`,
confirm that enough metacells remain for Pando, and only then scale workers.
`upstream_workers` parallelizes strata and `layer2_workers` parallelizes
metacell scoring; worker counts affect runtime and memory, not model semantics.

## Explicit setup

Use the same main workflow when the GEM or medium must be inspected or
customized first:

```r
gem <- rc_prepare_human2_gem(version = "2.0.0")
medium <- rc_make_medium_scenarios(
  gem,
  scenario = "compass_model_bounds"
)

result <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  medium_scenarios = medium
)
```

Keep one fixed metacell `gamma` across strata. Strata below the minimum
metacell count are recorded as skipped and excluded from calibration and
scoring. Advanced settings remain available through `metacell_args`,
`pando_args`, `layer1_args`, and `layer2_args`; defaults define the supported
main path.

The primary outputs are `result$layer1`, `result$grn_meta_modules`, and
`result$microcompass`. See the [workflow](docs/workflow.md) and
[public functions](docs/functions.md) for the compact contract.
