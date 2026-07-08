# RegCompassR strict multiome workflow

RegCompassR 的正式 multiome workflow 固定为：

```text
Seurat RNA+ATAC object
→ sample_id × condition × cell_type SuperCell metacells
→ metacell fragment aggregation
→ cell-type-stratified metacell LinkPeaks
→ RNA-GPR C_rel + ATAC-supported reaction_confidence
→ strict cached microCOMPASS
```

正式 workflow 不接受用户提供的 `peak_gene_links`，不复用 single-cell links，也不会在 fragment aggregation 失败后继续宣称 multiome analysis。

## Minimal workflow

```r
library(RegCompassR)

# 1. Validate a Seurat/Signac RNA+ATAC object.
rc_validate_multiome_input(
  object,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type"
)

# 2. Prepare a GEM. If save_rds is absent, a pinned Human-GEM YAML is downloaded,
# converted, validated, and cached.
gem <- rc_prepare_human2_gem(
  version = "2.0.0",
  save_rds = "Human2_2.0.0_regcompass.rds"
)

# Optional: add curated reaction roles and medium scenarios for microCOMPASS.
gem <- rc_annotate_reaction_roles(gem, reaction_role_table = reaction_roles)
medium <- rc_make_medium_scenarios(gem, scenario = "blood_like")

# 3. Run formal metacell-level multiome Layer 1.
# fragment_files must be valid fragments.tsv.gz files with indexes.
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
  min_metacells_for_linkpeaks = 80,
  and_method = "boltzmann",
  tau = 0.20
)

# 4. Select target reactions from Layer 1 evidence.
targets <- rc_select_target_reactions(
  layer1,
  method = "balanced_top_capacity",
  selection_mode = "balanced_rank",
  group_cols = c("condition", "cell_type"),
  top_n = 100,
  min_C_rel = 0.15,
  min_confidence = 0.25
)

# 5. Run strict cached microCOMPASS.
res <- rc_run_microcompass(
  layer1 = layer1,
  gem = gem,
  target_reactions = targets$reaction_id,
  medium_scenarios = medium,
  unit = "sample_celltype",
  target_direction = "both",
  omega = 0.95,
  solver = "highs",
  time_limit = 60
)

# 6. Optional differential testing and export.
stat <- rc_test_microcompass_differential(
  res,
  formula = score ~ condition,
  method = "lm",
  min_samples_per_group = 3
)

rc_export_microcompass(res, "RegCompassR_run")
```

## Function guide

| Step | Function | Notes |
|---|---|---|
| Input check | `rc_validate_multiome_input()` | Checks Seurat/Signac assays and required metadata. |
| GEM | `rc_prepare_human2_gem()` | Loads cached RDS or downloads/converts pinned Human-GEM YAML. |
| Roles | `rc_annotate_reaction_roles()` | Optional, improves micro-GEM role annotations. |
| Medium | `rc_make_medium_scenarios()` | Optional named medium scenario table. |
| Formal Layer 1 | `rc_run_regcompass_multiome_metacell()` | Builds sample-aware metacells, requires fragment aggregation, recomputes metacell LinkPeaks by `link_stratum_cols`, and returns Layer 1 evidence. |
| Target selection | `rc_select_target_reactions()` | Selects reactions using `C_rel` and reaction confidence. |
| microCOMPASS | `rc_run_microcompass()` | Runs strict cached target-local LP scoring. |
| Statistics | `rc_test_microcompass_differential()` | Optional sample/group-level testing. |
| Export | `rc_export_microcompass()` | Writes strict microCOMPASS outputs only. |

## Layer 1 outputs

`rc_run_regcompass_multiome_metacell()` returns a Layer 1 object with:

- `C_rel`
- `reaction_confidence`
- `metacell_meta` / `pool_meta`
- `rna_metacell_logcpm`
- `rna_metacell_detection`
- `metacell_peak_gene_links`
- `peak_gene_link_source`

## Exported files

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
