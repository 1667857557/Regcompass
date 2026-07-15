# RegCompassR

RegCompassR 1.4 implements a two-stage, shared-structure multiome workflow for directional reaction-potential analysis.

## Workflow

```text
Seurat RNA+ATAC object
→ define retained condition × sample × cell-type strata
→ upstream parallel stage: one worker per stratum
     Meta cell
     → fragment aggregation
     → metacell LinkPeaks
     → Layer 1 RNA/ATAC-GPR evidence
     → Pando GRN
     → metabolic GRN projection
     → meta-module expansion
→ global barrier: every retained stratum must complete
→ merge all metacells and rerun expression-capacity calibration globally
→ union all stratum meta-modules
→ build one shared FASTCORE-completed GEM per medium scenario
→ release upstream workers
→ downstream parallel stage: shared target-direction × metacell LP tasks
```

The downstream model has the same stoichiometric matrix, bounds, target definitions, and medium for every metacell. Each metacell retains its own evidence-derived penalty vector.

## Installation

```r
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

## Integrated analysis

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

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

result <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "RegCompassR_v1.4",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  model_mode = "meta_module_gem",
  medium_scenarios = medium,
  metacell_args = list(
    gamma = 150,
    adaptive_gamma = TRUE,
    min_cells_pre_metacell = 100,
    min_metacell_size = 10,
    min_metacells_post_metacell = 10
  ),
  pando_args = list(
    min_metacells = 20,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr",
      parallel = FALSE
    ),
    padj_threshold = 0.05,
    min_model_rsq = 0.1,
    top_k_neighbors = 5,
    expansion_mode = "ordered_once"
  ),
  calibration_args = list(
    bootstrap = FALSE,
    promiscuity_mode = "sqrt",
    and_method = "boltzmann",
    tau = 0.20,
    or_method = "sum_sqrtK"
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    parallel = TRUE,
    model_params = list(
      strict = TRUE,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000
    )
  ),
  BPPARAM_upstream = rc_default_bpparam(workers = 6, backend = "snow"),
  BPPARAM_layer2 = rc_default_bpparam(workers = 12, backend = "snow")
)
```

`condition`-specific medium bounds are rejected because they would make the structural model differ between conditions. Use `condition = "all"` when cross-condition scores must share one GEM.

## Barrier and outputs

No global calibration, union GEM construction, or LP scoring starts unless every retained stratum completes all upstream steps. Diagnostics are written to:

```text
upstream_expected_strata.tsv.gz
upstream_stratum_status.tsv.gz
```

Main outputs:

```r
result$layer1$C_raw
result$layer1$C_rel
result$grn_meta_modules$core_gene_reaction
result$grn_meta_modules$reaction_membership
result$microcompass$score
result$microcompass$penalty
result$microcompass$vmax
result$microcompass$model_cache_summary
result$microcompass$lp_diagnostics
```

Score rows use a shared identity without sample- or module-specific model labels:

```text
reaction=<reaction>::direction=<direction>::medium=<medium>::condition=all
```

## Public API

`rc_run_regcompass()` is the supported end-to-end workflow. Former public stage-specific runners and cache constructors are now internal implementation details to avoid competing orchestration paths.

## Interpretation

Reaction capacities and microCOMPASS scores are model-based reaction potentials, not measured fluxes. Metacells are technical aggregation units and are not biological replicates. Results remain conditional on the GEM, GPR mapping, reaction-role annotations, medium constraints, GRN inference, and penalty model.
