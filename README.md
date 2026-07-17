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
  celltype_col = "cell_type"
)
```

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
