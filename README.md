# RegCompassR

RegCompassR provides one supported RNA+ATAC workflow:

```text
condition × sample × cell type
→ metacells and stratum-specific Pando GRNs
→ local FASTCORE meta-modules
→ sample-balanced Q95 diagnostics and one shared GEM
→ directional microCOMPASS scoring
```

## Installation

```r
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

Fragment-enabled runs additionally require a MACS2/MACS3 executable. Pass its
path through `metacell_args$macs2_path` when it is not available as `macs2` on
`PATH`.

## Quick start

The one-shot entry point prepares Human-GEM 2 and the shared model-bound medium
when they are not supplied. This is the default human path; set
`species = "mouse"` to route setup to Mouse-GEM and the mouse physiological
medium:

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
  fragment_files = FALSE,  # use existing ATAC peak counts; pass paths to re-aggregate fragments
  species = "human",  # default; use "mouse" for Mouse-GEM + mouse medium
  gem_version = "2.0.0",
  medium_scenario = "physiologic",
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

Set `fragment_files = FALSE` when no matching fragment files are available; the
workflow skips fragment aggregation and carries the object's ATAC peak raw counts
into metacell and Pando analysis. When matching fragments are available, pass a
path, named list, or manifest in `fragment_files`; RegCompass re-aggregates ATAC
peak raw counts from fragments before downstream analysis. For mouse data, use
the matching genome and set `species = "mouse"`; the one-shot setup then
prepares Mouse-GEM 1.8.0 and `mouse_plasma` through the `"physiologic"` medium
shortcut. To use a fully custom medium table, build it with
`rc_make_medium_scenarios()` and pass it as `medium_scenarios`, which overrides
`medium_scenario`.

## Choosing analysis parameters

The values above are starting points, not fixed analysis parameters. Choose
them before a full run from stratum sizes, expected biological heterogeneity,
and a small pilot run.

| Setting | Role | Selection guidance |
| --- | --- | --- |
| `metacell_args$gamma` | Approximate cells represented by each metacell | Use a smaller value for more metacells and finer heterogeneity, or a larger value for stronger aggregation. Keep one value across all strata so resolution is comparable. |
| `metacell_args$min_cells_per_stratum` | Minimum cells in each condition × sample × cell-type stratum | Set high enough to avoid unstable strata, but check that every biological sample retains at least one stratum. |
| `metacell_args$min_metacell_size` | Flags undersized, low-power metacells | Increase when sparse RNA/ATAC profiles are unstable; do not use it to compensate for an unsuitable `gamma`. |
| `metacell_args$macs2_path` | MACS2/MACS3 executable for fragment-enabled runs | Required when the executable is not discoverable as `macs2` on `PATH`. |
| `metacell_args$peak_calling_effective_genome_size` | MACS effective genome size | Use a species-matched value; defaults are inferred from annotated `hg*`/`GRCh*` or `mm*`/`GRCm*` peak ranges. |
| `metacell_args$peak_calling_args` | Additional `Signac::CallPeaks()` arguments | Keep one policy across strata. Use only justified MACS options and record them with the run. |
| `pando_args$min_metacells` | Minimum metacells required for Pando | It must be compatible with `floor(n_cells / gamma)`. Strata below it are skipped, so inspect `00_strata/stratum_workflow_status.tsv.gz` after a pilot. |
| `pando_infer_args` | Pando model and correlation/FDR filters | Start with the shown GLM settings; tighten correlation thresholds only when the pilot produces excessive weak edges. Apply one policy to every stratum. |
| `layer1_args$local_fastcore` | Completes each local metabolic module | Keep enabled for the canonical path. |
| `layer1_args$sample_balance` | Defines the sampling estimand for Q95 and relative-state diagnostics | Keep `TRUE` for biological-replicate inference. Every sample receives equal total mass globally, and weights are recomputed inside each Q95 stratum. Absolute metacell activity is not rescaled. |
| `layer1_args$expression_batch_correction` | Optional technical-batch correction | Keep `"none"` unless a documented technical batch exists. If using `"limma"`, provide technical and preserved biological design columns; never remove `sample_id` as batch. |
| `layer2_args$target_direction` | Forward, reverse, or both-direction scoring | Use `"both"` unless the GEM direction or scientific question justifies one direction. |
| `layer2_args$solver` | LP solver | `"highs"` is the default open-source choice; use Gurobi only in a licensed environment. |

`sample_balance = TRUE` uses \(w_{si}=1/(S n_s)\) for global diagnostics and
recomputes \(w_{sci}=1/(S_c n_{sc})\) inside each Q95 stratum. These are
sampling weights, not enzyme-activity multipliers. The primary absolute
reaction support remains the zero-preserving metacell value, while the default
`sample_celltype` Layer 2 unit gives each biological sample one inference unit
per represented cell type.

A practical pilot is to tabulate cells per strict stratum, choose one `gamma`,
confirm that enough metacells remain for Pando, and only then scale workers.
`upstream_workers` parallelizes strata and `layer2_workers` parallelizes
metacell scoring; worker counts affect runtime and memory, not model semantics.

## Explicit setup

Use the same main workflow when the GEM or medium must be inspected or
customized first:

```r
human2_gem <- rc_prepare_human2_gem(version = "2.0.0")
# Or choose Mouse-GEM explicitly when analyzing mouse data:
# mouse_gem <- rc_prepare_mouse_gem(version = "1.8.0")
medium <- rc_make_medium_scenarios(
  human2_gem,
  scenario = "compass_model_bounds"
)

result <- rc_run_regcompass(
  object = object,
  gem = human2_gem,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  species = "human",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  medium_scenarios = medium
)
```

## Published human medium presets

`compass_model_bounds` remains the technical default. Human-only biological
backgrounds are selected explicitly:

```r
plasma <- rc_make_medium_scenarios(
  gem,
  scenario = "normal_human_plasma"
)

high_glucose <- rc_make_medium_scenarios(gem, scenario = "high_glucose")
low_glucose <- rc_make_medium_scenarios(gem, scenario = "low_glucose")
high_lactate <- rc_make_medium_scenarios(gem, scenario = "high_lactate")
low_lactate <- rc_make_medium_scenarios(gem, scenario = "low_lactate")
rpmi <- rc_make_medium_scenarios(gem, scenario = "rpmi1640")
```

The presets close uptake for exchanges not represented by the selected
background and reopen listed nutrients. They do not convert concentration in
mM into physical flux. Glucose and lactate concentrations define relative
sensitivity caps; other listed metabolites are treated as available. Every
preset row records the human paper citation, DOI, PMID, concentration and bound
provenance. See `?rc_make_medium_scenarios` and
[`docs/functions.md`](docs/functions.md) for the exact values and references.

Users can supply either exact reaction bounds with `custom_medium` or a
metabolite availability table with `custom_metabolites`:

```r
custom <- rc_make_medium_scenarios(
  gem,
  scenario = "custom",
  custom_metabolites = data.frame(
    metabolite_name = c("glucose", "lactate"),
    metabolite_pattern = c("glucose|glc", "lactate|lactic acid"),
    available = TRUE,
    concentration_mM = c(3, 8),
    uptake_fraction = c(0.12, 0.40),
    target_exchange_flag = TRUE,
    required_match = TRUE,
    reference_doi = "project-specific reference"
  )
)
```

Keep one fixed metacell `gamma` across strata. Strata below the minimum
metacell count are recorded as skipped and excluded from calibration and
scoring. Advanced settings remain available through `metacell_args`,
`pando_args`, `layer1_args`, and `layer2_args`; defaults define the supported
main path.

The primary outputs are `result$layer1`, `result$grn_meta_modules`, and
`result$microcompass`. See the [workflow](docs/workflow.md) and
[public functions](docs/functions.md) for the compact contract. The
[`regcompass-workflow` vignette](vignettes/regcompass-workflow.Rmd) provides a
complete, input-to-output walkthrough with explicit and one-shot entry points.
