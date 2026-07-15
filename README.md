# RegCompassR

RegCompassR converts RNA+ATAC metacells, GPR evidence, Pando regulatory networks, Human-GEM structure, and medium constraints into directional microCOMPASS reaction-potential scores.

## Current workflow

```text
Seurat RNA+ATAC object
→ split by condition × sample × cell type
→ one upstream worker per strict stratum:
     SuperCell2 metacells
     fragment aggregation
     Pando internal peak-gene and TF-gene modeling
     Pando-derived metacell ATAC confidence
     Pando GRN
     GRN-to-reaction mapping and meta-module expansion
→ wait for every retained stratum and every sample to complete
→ validate the all-strata artifact barrier
→ stop and release the upstream worker pool
→ combine Human-GEM GPR-gene logCPM from all metacells
→ recompute one global gene-score scale and reaction-capacity matrix
→ recompute one reaction-wise Q95 scale across all metacells
→ union all stratum meta-modules
→ complete one shared global meta-module GEM with add-only FASTCORE
→ align all metacells to the same reaction universe
→ compute one penalty vector per metacell
→ start a fresh worker pool for shared-model/medium × metacell tasks
→ load the shared GEM once per metacell task and evaluate all target directions
```

No global recalibration or GEM construction occurs from a partial set of successful strata. If one retained stratum fails, the workflow writes `00_strata/upstream_barrier.tsv.gz`, releases the upstream workers, and stops.

The shared GEM is the structural reference for every metacell. Biological differences enter the LP objective through metacell-specific penalties, not through sample-specific stoichiometric models. After all upstream artifacts are complete, RegCompass combines GPR-gene logCPM across all metacells, recomputes gene scores and GPR reaction capacities on that common population, and then applies one reaction-wise Q95 calibration. Stratum-local `C_raw` and `C_rel` values are not used for cross-sample scoring; stratum-specific LinkPeaks evidence remains in the reaction-confidence term.

## Structural modes

| Mode | Structural model |
|---|---|
| `meta_module_gem` | Union of all strict-stratum Pando meta-modules, completed once with Human-GEM support reactions. Recommended. |
| `full_gem` | One complete medium-constrained Human-GEM shared by all metacells. |

Condition-specific medium bounds are rejected by the integrated shared-GEM workflow because they would create different feasible spaces across conditions.

## Installation

```r
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

The installed Pando package must retain repository metadata pointing to `1667857557/Pando_regcompass`.

## Prepare Human-GEM and medium constraints

```r
library(RegCompassR)

gem <- rc_prepare_human2_gem_v12(version = "2.0.0")

medium <- data.frame(
  medium_scenario_id = "brain_tumor_medium",
  exchange_reaction_id = curated_exchange_ids,
  lb = curated_uptake_lower_bounds,
  ub = rep(1000, length(curated_exchange_ids)),
  available = TRUE,
  condition = "all",
  stringsAsFactors = FALSE
)
```

Reaction directions are derived from active numerical bounds:

```text
lb < 0 < ub     reversible
0 <= lb < ub    forward-only
lb < ub <= 0    reverse-only
lb = ub = 0     blocked
```

## Integrated analysis

```r
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

result <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "RegCompass_global_metacell",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  rna_assay = "RNA",
  atac_assay = "ATAC",
  model_mode = "meta_module_gem",
  medium_scenarios = medium,
  metacell_args = list(
    gamma = 150,
    adaptive_gamma = TRUE,
    min_cells_per_stratum = 100,
    min_metacell_size = 10,
    min_metacells_per_stratum = 10
  ),
  layer1_args = list(
    promiscuity_mode = "sqrt",
    and_method = "boltzmann",
    tau = 0.20
  ),
  pando_args = list(
    min_metacells = 20,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr"
    ),
    padj_threshold = 0.05,
    min_model_rsq = 0.1,
    top_k_neighbors = 5,
    min_shared_tfs = 1,
    expansion_mode = "ordered_once"
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      strict = TRUE,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000
    )
  ),
  upstream_workers = 6,
  layer2_workers = 12,
  parallel_backend = "snow"
)
```

The integrated runner fixes `unit = "metacell"`, uses one shared structural GEM, and prevents `layer2_args` from overriding those invariants.

## Main outputs

```r
result$upstream_status
result$upstream_barrier
result$layer1$rna_metacell_logcpm
result$layer1$global_gene_score
result$layer1$C_raw
result$layer1$C_rel
result$layer1$capacity_calibration_scope
result$layer1$reaction_confidence
result$grn_meta_modules$core_gene_reaction
result$grn_meta_modules$reaction_membership
result$grn_meta_modules$global_core_reactions
result$grn_meta_modules$global_reaction_membership
result$microcompass$score
result$microcompass$penalty
result$microcompass$vmax
result$microcompass$feasible
result$microcompass$model_cache_summary
result$microcompass$lp_diagnostics
```

The complete result is saved as:

```text
RegCompass_global_metacell/regcompass_global_metacell_result.rds
```

## Retired public APIs

The former staged integrated entry point and sample-specific meta-module cache are no longer exported:

```text
rc_run_regcompass_multiome_metacell
rc_build_meta_module_gem_cache
rc_load_metacell_object_from_run
```

Use `rc_run_regcompass()` for the supported workflow. Lower-level metacell, Pando, FASTCORE, and LP functions remain available where they are still part of the current architecture.

## Interpretation limits

- FASTCORE is an LP-based compact reconstruction algorithm, not an exact minimum-cardinality MILP.
- Steady-state feasibility does not establish thermodynamic or kinetic feasibility.
- Reaction capacity and microCOMPASS scores are multiome-supported reaction potentials, not measured fluxes.
- Metacells are technical aggregation units and are not independent biological replicates.
- Cross-condition inference must aggregate or model metacells at the biological-sample level.
