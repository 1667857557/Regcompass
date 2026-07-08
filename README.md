# RegCompassR workflow

RegCompassR runs a strict multiome workflow:

```text
Seurat RNA+ATAC object
→ sample-aware metacells and fragments
→ metacell LinkPeaks
→ Layer 1 reaction capacity/confidence
→ strict microCOMPASS
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
  gamma = 100,
  min_metacell_size = 20,
  link_stratum_cols = "cell_type",
  min_metacells_for_linkpeaks = 10,
  future_plan = "sequential"
)

targets <- rc_select_target_reactions(
  layer1,
  method = "balanced_top_capacity",
  selection_mode = "balanced_rank",
  group_cols = c("condition", "cell_type"),
  top_n = 500,
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
  solver = "highs"
)

stat <- rc_test_microcompass_differential(
  res,
  formula = score ~ condition,
  method = "lm"
)

rc_export_microcompass(res, "RegCompassR_run")
```

`fragment_files` may be a single fragment path or a named vector/list such as `list(ATAC = "fragments.tsv.gz")`.

## Function checklist

| Function | Required role in the example |
|---|---|
| `rc_validate_multiome_input()` | Check Seurat RNA/ATAC assays and metadata columns. |
| `rc_prepare_human2_gem()` | Load a pinned Human2 GEM with `gpr_table`. |
| `rc_annotate_reaction_roles()` | Add reaction roles used by medium and microCOMPASS steps. |
| `rc_make_medium_scenarios()` | Build the medium table passed as `medium_scenarios`. |
| `rc_run_regcompass_multiome_metacell()` | Build metacells, aggregate fragments, recompute LinkPeaks, and return Layer 1. |
| `rc_select_target_reactions()` | Select Layer 2 target reactions from Layer 1 results. |
| `rc_run_microcompass()` | Run target-local LP scoring. |
| `rc_test_microcompass_differential()` | Test score differences using `result$unit_meta`. |
| `rc_export_microcompass()` | Write matrices, diagnostics, and `session_info.txt`. |

## Main outputs

- `layer1`: includes `C_rel`, `reaction_confidence`, `metacell_meta`, `pool_meta`, `rna_metacell_logcpm`, `rna_metacell_detection`, and `metacell_peak_gene_links`.
- `res`: includes `score`, `penalty`, `vmax`, `feasible`, `medium_scenarios`, diagnostics, and `unit_meta`.
- `stat`: differential-test table from `rc_test_microcompass_differential()`.
