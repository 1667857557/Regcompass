# RegCompassR workflow

RegCompassR runs a strict multiome workflow:

```text
Seurat RNA+ATAC object
→ pre-filter condition × sample × cell type strata
→ build SuperCell2 RNA+ATAC metacells inside each retained stratum
→ post-filter strata by actual metacell count
→ aggregate ATAC fragments separately for each retained stratum
→ recompute stratum-specific metacell LinkPeaks
→ Layer 1 reaction capacity/confidence
→ Layer 2 GEM scoring
```

## Minimal example

```r
library(RegCompassR)

rc_validate_multiome_input(
  object,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type"
)

gem <- rc_prepare_human2_gem(version = "2.0.0")
gem <- rc_annotate_reaction_roles(gem, reaction_role_table = reaction_roles)
medium <- rc_make_medium_scenarios(gem, scenario = "blood_like")

bp <- rc_default_bpparam(workers = 4, backend = "snow")

layer1 <- rc_run_regcompass_multiome_metacell(
  object = object,
  gpr_table = gem$gpr_table,
  outdir = "RegCompassR_run",
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  rna_assay = "RNA",
  atac_assay = "ATAC",
  gamma = 150,
  adaptive_gamma = TRUE,
  min_cells_pre_metacell = 100,
  min_metacell_size = 10,
  min_metacells_post_metacell = 10,
  future_plan = "sequential",
  BPPARAM_metacell = FALSE,      # keep SuperCell2/Seurat object work serial by default
  BPPARAM_linkpeaks = bp,        # LinkPeaks strata can be parallelized
  BPPARAM_layer1 = bp            # Layer 1 reaction/bootstrap loops can be parallelized
)

targets <- rc_select_target_reactions(
  layer1,
  method = "balanced_top_capacity",
  selection_mode = "balanced_rank",
  group_cols = c("condition", "cell_type"),
  top_n = 100,
  min_C_rel = 0.15,
  min_confidence = 0.25
)

res <- rc_run_microcompass(
  layer1 = layer1,
  gem = gem,
  target_reactions = targets$reaction_id,
  medium_scenarios = medium,
  unit = "sample_celltype",
  target_direction = "both",
  # Default: COMPASS-style full GEM per medium scenario.
  solver = "highs",
  parallel = TRUE,
  BPPARAM = bp
)

stat <- rc_test_microcompass_differential(
  res,
  formula = score ~ condition,
  method = "lm"
)

rc_export_microcompass(res, "RegCompassR_run")
```

`fragment_files` may be a single fragment path for one-sample data or a manifest with `sample_id`, `assay`, and `fragment_file` columns for multi-sample data. The formal workflow now aggregates fragments after post-metacell filtering and writes one set of metacell fragment outputs per retained `condition × sample × cell type` stratum.

## Function checklist

| Function | Required role in the example |
|---|---|
| `rc_validate_multiome_input()` | Check Seurat RNA/ATAC assays and metadata columns. |
| `rc_prepare_human2_gem()` | Load a pinned Human2 GEM with `gpr_table`. |
| `rc_annotate_reaction_roles()` | Add reaction roles used by medium and microCOMPASS steps. |
| `rc_make_medium_scenarios()` | Build the medium table passed as `medium_scenarios`. |
| `rc_run_regcompass_multiome_metacell()` | Run the current formal workflow: strict-stratum metacells, post-filtered stratum-wise fragment aggregation, parallelizable metacell LinkPeaks, and Layer 1. It can also run Layer 2 when fixed targets are supplied through `layer2_gem` and `layer2_target_reactions`. |
| `rc_select_target_reactions()` | Select Layer 2 target reactions from Layer 1 results. |
| `rc_run_microcompass()` | Run COMPASS-style LP scoring; by default it caches one full GEM per medium scenario and maps selected target reactions to that full-model input. Use `parallel = TRUE` with `BPPARAM` for GEM task parallelism. |
| `rc_build_full_gem_cache()` | Optional lower-level helper used by the default `rc_run_microcompass()` strategy (`microgem_params = list(strategy = "full_gem")`). |
| `rc_build_module_gem_cache()` | Optional lower-level helper for the older module-level strategy (`microgem_params = list(strategy = "module_meso_gem")`). |
| `rc_test_microcompass_differential()` | Test score differences using `result$unit_meta`. |
| `rc_export_microcompass()` | Write matrices, diagnostics, and `session_info.txt`. |

## Strict metacell strata and LinkPeaks gates

The formal metacell workflow uses one fixed analysis stratum definition for
metacell construction, post-metacell filtering, fragment aggregation, and LinkPeaks:

```text
condition × sample × cell type
```

Input validation only confirms that required assays, reductions, barcodes, and
metadata columns are usable. It does not guarantee every stratum has enough data
for downstream analysis. The workflow applies two auditable hard filters:

1. Before SuperCell, each `condition × sample × cell type` stratum must have at
   least `min_cells_pre_metacell = 100` original cells. Strata below this
   threshold are excluded before metacell construction.
2. After metacell construction, each same stratum must have at least
   `min_metacells_post_metacell = 10` actual generated metacells. Strata below
   this threshold are removed from all downstream bundles.
3. After post-metacell filtering, ATAC fragments are aggregated independently for
   each retained strict stratum. These metacell fragment files are then registered
   on the merged metacell Seurat object for Signac.
4. LinkPeaks is recomputed independently within the same strict stratum and uses
   the same `min_metacells_post_metacell` threshold; any retained stratum with
   fewer actual metacells is treated as an internal invariant failure. This
   threshold is independent of Signac's `min.cells` detection parameter. In
   Signac, `min.cells` is a feature-detection gate for peaks/genes observed in
   at least that many cells/metacells, so RegCompassR passes `min.cells = 3` by
   default during metacell LinkPeaks to avoid dropping sparsely detected genes
   solely because they are present in fewer than 10 metacells.
5. Strata excluded by either gate do not enter fragment aggregation, LinkPeaks,
   Layer 1, or microCOMPASS. Excluded cells/metacells are retained only in QC
   reports under `00_stratum_qc/`.

`state_col` can be kept in metadata for later summaries, but it does not change
the formal metacell/LinkPeaks stratum in this workflow.

## Current API notes

Use the formal entry point `rc_run_regcompass_multiome_metacell()` for new analyses. The older convenience wrappers `rc_make_metacells()` and `rc_import_metacells()` have been removed to avoid conflicting filtering semantics; use `rc_make_supercell2_metacells()` / `rc_import_supercell2_metacells()` only for lower-level debugging or custom workflows. Deprecated aliases such as `min_cells_per_stratum`, `min_metacells_per_stratum`, `min_metacells_for_linkpeaks` on the formal entry point, and `filter_low_power_metacells` are no longer part of the current API.

## Layer 2 GEM cache strategies

`rc_run_microcompass()` now defaults to a COMPASS-style full-GEM analysis. With
the default `microgem_params = list()` setting, RegCompassR builds one validated
full GEM per `medium_scenario_id`, applies the matching medium constraints, and
uses that complete model when scoring each selected target reaction. This avoids
selecting reactions by metabolic module unless you explicitly request it.

Use these `microgem_params$strategy` values when you need a different structural
cache:

| Strategy | Behavior |
|---|---|
| `"full_gem"` | Default. Cache the complete GEM per medium scenario, matching standard COMPASS-style full-model LP inputs. |
| `"module_meso_gem"` | Cache one module-level meso-GEM per module and medium scenario; requires `gem$reaction_meta$metabolic_module` or a custom `module_col`. |
| `"target_khop"` | Build target-local k-hop micro-GEMs through `rc_build_microgem_cache()`. |
| `"auto"` | Try target-local k-hop micro-GEMs and fall back to module meso-GEMs when strict closure diagnostics fail. |

For example, to reproduce the previous module-based workflow:

```r
res_module <- rc_run_microcompass(
  layer1 = layer1,
  gem = gem,
  target_reactions = targets$reaction_id,
  medium_scenarios = medium,
  unit = "sample_celltype",
  target_direction = "both",
  microgem_params = list(strategy = "module_meso_gem"),
  solver = "highs"
)
```

## Main outputs

- `layer1`: includes `C_rel`, `reaction_confidence`, `metacell_meta`, `unit_meta`, `rna_metacell_logcpm`, `rna_metacell_detection`, and `metacell_peak_gene_links`.
- `res`: includes `score`, `penalty`, `vmax`, `feasible`, `medium_scenarios`, diagnostics, and `unit_meta`.
- `stat`: differential-test table from `rc_test_microcompass_differential()`.
