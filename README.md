# RegCompassR workflow

```text
Seurat RNA+ATAC object
→ sample-aware metacells
→ metacell fragment aggregation
→ metacell LinkPeaks
→ RNA-GPR capacity + ATAC-supported confidence
→ strict cached microCOMPASS
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
  link_stratum_cols = "cell_type",
  min_metacells_for_linkpeaks = 20
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
  solver = "highs"
)

stat <- rc_test_microcompass_differential(
  res,
  formula = score ~ condition,
  method = "lm"
)

rc_export_microcompass(res, "RegCompassR_run")
```

## Function checklist

| Function | Use |
|---|---|
| `rc_validate_multiome_input()` | Validate Seurat assays and metadata. |
| `rc_prepare_human2_gem()` | Load or build a pinned Human2 GEM with `gpr_table`. |
| `rc_annotate_reaction_roles()` | Add curated/inferred reaction roles. |
| `rc_make_medium_scenarios()` | Create medium scenario table for microCOMPASS. |
| `rc_run_regcompass_multiome_metacell()` | Build metacells, aggregate fragments, recompute LinkPeaks, and return Layer 1 evidence. |
| `rc_select_target_reactions()` | Choose target reactions from `C_rel` and reaction confidence. |
| `rc_run_microcompass()` | Run strict cached target-local LP scoring. |
| `rc_test_microcompass_differential()` | Test score differences by sample metadata. |
| `rc_export_microcompass()` | Write strict microCOMPASS matrices and diagnostics. |

## Key outputs

`layer1` includes `C_rel`, `reaction_confidence`, `metacell_meta`, `pool_meta`, `rna_metacell_logcpm`, `rna_metacell_detection`, `metacell_peak_gene_links`, and `peak_gene_link_source`.

`res` includes `score`, `penalty`, `vmax`, `feasible`, `target_direction`, `medium_scenarios`, `microgem_cache_summary`, `microgem_diagnostics`, `lp_diagnostics`, `penalty_components`, and `unit_meta`.

## Export layout

```text
02_medium/medium_scenarios.tsv.gz
02_medium/medium_sensitivity_summary.tsv.gz
03_microgem/closure_diagnostics.tsv.gz
03_microgem/microgem_cache_summary.tsv.gz
04_microcompass/strict_score_matrix.rds
04_microcompass/strict_penalty_matrix.rds
04_microcompass/vmax_matrix.rds
04_microcompass/feasible_matrix.rds
04_microcompass/penalty_components.rds
04_microcompass/lp_diagnostics.tsv.gz
session_info.txt
```
