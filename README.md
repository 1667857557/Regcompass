# RegCompassR microCOMPASS quick guide

microCOMPASS builds structural target-local micro-GEMs from GEM topology, target reactions, reaction roles, and medium scenarios. RNA+ATAC-GPR evidence is used only for unit-specific LP penalties; it does not prune reactions.

Output means multiome-supported reaction capacity potential. It is not true flux, enzyme activity, uptake/secretion flux, ATAC causality, or in vivo medium truth.

## Minimal workflow

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

mc <- rc_make_metacells(
  object = object,
  outdir = "RegCompassR_run/00_metacells",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  rna_assay = "RNA",
  atac_assay = "ATAC",
  gamma = 75,
  min_cells_per_stratum = 100,
  min_metacell_size = 20,
  filter_low_power_metacells = TRUE
)

gem <- rc_prepare_human2_gem(
  version = "2.0.0",
  save_rds = "Human2_2.0.0_regcompass.rds"
)
gem <- rc_annotate_reaction_roles(gem, reaction_role_table = reaction_roles)

medium <- rc_make_medium_scenarios(
  gem,
  scenario = c("blood_like", "low_glucose", "low_glutamine", "lactate_available")
)

layer1 <- rc_run_layer1_multiome(
  gpr_table = gem$gpr_table,
  rna_metacell_counts = mc$rna_counts,
  metacell_meta = mc$metacell_meta_used,
  atac_metacell_counts = mc$atac_counts,
  peak_gene_links = peak_gene_links,
  stratum_col = "cell_type",
  and_method = "boltzmann",
  tau = 0.20,
  q = 0.95,
  q95_n0 = 80
)

targets <- rc_select_target_reactions(
  layer1,
  method = "balanced_top_capacity",
  selection_mode = "balanced_rank",
  group_cols = c("condition", "cell_type"),
  top_n = 100,
  top_n_per_group = 30,
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
  omega = 0.95,
  use_gapfilled_for_score = FALSE,
  parallel = TRUE,
  solver = "highs",
  time_limit = 60
)

stat <- rc_test_microcompass_differential(
  res,
  formula = score ~ condition,
  method = "lm",
  min_samples_per_group = 3
)

rc_export_microcompass(res, "RegCompassR_run")
```

## Function checklist

| Function | Purpose | Returns |
|---|---|---|
| `rc_validate_multiome_input()` | Validate Seurat/Signac assays and metadata. | Invisible `TRUE`; errors on invalid input. |
| `rc_make_metacells()` | Build sample/condition/cell-type-aware RNA+ATAC metacells. | Counts, metadata, membership, diagnostics. |
| `rc_import_metacells()` | Import saved metacell outputs. | Same structure as `rc_make_metacells()`. |
| `rc_prepare_human2_gem()` | Load pinned preconverted Human2 RDS; no automatic conversion. | Validated Human2 GEM. |
| `rc_read_gem()` | Read generic GEM RDS with `model_info` by default. | Validated GEM. |
| `rc_annotate_reaction_roles()` | Add curated/inferred reaction roles. | GEM with `reaction_roles`. |
| `rc_make_medium_scenarios()` | Build named exchange-bound scenarios. | Medium scenario table. |
| `rc_run_layer1_multiome()` | Compute RNA+ATAC-GPR reaction evidence. | `C_rel`, confidence, diagnostics, unit metadata. |
| `rc_select_target_reactions()` | Select target reactions only. | Target table with balanced ranking diagnostics. |
| `rc_build_target_microgem()` | Build one structural target-local micro-GEM. | Micro-GEM with closure/medium diagnostics. |
| `rc_build_microgem_cache()` | Cache micro-GEMs by target direction and medium scenario. | Named cache plus `microgem_cache_summary`. |
| `rc_compute_multiome_penalty()` | Convert unit evidence to LP penalties. | Penalty matrix, components, `evidence_policy = "penalty_only"`. |
| `rc_run_microcompass()` | Run cached strict two-step LP per target/scenario/unit. | Score, penalty, vmax, feasibility, LP diagnostics. |
| `rc_test_microcompass_differential()` | Test sample-level scores with `lm`, `wilcox`, or `limma_continuous`. | Differential result table. |
| `rc_export_microcompass()` | Write standardized outputs. | Files under `02_medium`, `03_microgem`, `04_microcompass`. |

## Required tables

### Human2 RDS

`rc_prepare_human2_gem()` expects a pinned preconverted RDS with `S`, `lb`, `ub`, `reaction_meta`, `metabolite_meta`, `gpr_table`, and `model_info`.

### Medium scenario table

Usually create this with `rc_make_medium_scenarios()`. A custom table needs at least:

```text
medium_scenario_id  exchange_reaction_id  lb   ub    available
blood_like          EX_glc_D_e            -10  1000  TRUE
low_glucose         EX_glc_D_e             -5  1000  TRUE
```

Optional columns include `metabolite_id`, `condition`, `evidence_source`, and `assumption_level`.

### Reaction role table

```text
reaction_id  role       role_source
EX_glc_D_e   exchange   curated
R_HEX1       internal   curated
```

Common roles include `internal`, `boundary_like`, `exchange`, `transport`, `demand`, `sink`, `cofactor_recycle`, and `unknown`.

## Exported outputs

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

## Interpretation

- Strict feasible + low penalty + stable medium sensitivity: stronger capacity-potential evidence.
- Strict infeasible: score is 0 for that structural micro-GEM and medium scenario.
- ATAC confidence means chromatin-expression concordance, not regulatory causality.
- DNA methylation is not modeled.
