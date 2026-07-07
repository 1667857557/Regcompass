# RegCompassR microCOMPASS quick guide

RegCompassR microCOMPASS uses structural target-local micro-GEMs defined by GEM topology, target reactions, reaction roles, and medium scenarios. RNA+ATAC-GPR evidence is used only to compute unit-specific penalties and does not remove reactions from the structural micro-GEM.

RegCompassR reports multiome-supported reaction capacity potential, not true metabolic flux, enzyme activity, uptake/secretion flux, ATAC causality, or in vivo medium truth.

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

gem <- rc_read_gem("Human2_2.0.0_regcompass.rds")
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

| Function | Use | Main output |
|---|---|---|
| `rc_validate_multiome_input()` | Check Seurat/Signac assays and metadata. | Invisible `TRUE`; stops on invalid input. |
| `rc_make_metacells()` | Build sample/condition/cell-type-aware RNA+ATAC metacells. | Counts, metadata, membership, diagnostics. |
| `rc_import_metacells()` | Import saved metacell directories. | Same structure as `rc_make_metacells()`. |
| `rc_read_gem()` | Read a GEM RDS; requires `model_info` by default. | Validated GEM. |
| `rc_prepare_human2_gem()` | Load a pinned preconverted Human2 RDS with provenance; automatic conversion is not bundled. | Validated Human2 GEM, or an error if no RDS/provenance is supplied. |
| `rc_annotate_reaction_roles()` | Add curated or inferred reaction roles. | GEM with `reaction_roles`. |
| `rc_make_medium_scenarios()` | Create named medium-scenario bounds. | Scenario table for exchanges. |
| `rc_run_layer1_multiome()` | Convert metacell RNA/ATAC and GPRs into reaction evidence. | `C_rel`, confidence, diagnostics, unit metadata. |
| `rc_select_target_reactions()` | Select targets without expanding networks. | Balanced ranked target table. |
| `rc_extract_cofactor_modules()` | Find cofactor-related candidate reactions from GEM annotations. | Cofactor module diagnostic table. |
| `rc_build_microgem_cache()` | Build structural micro-GEMs once per target direction and medium scenario. | Named micro-GEM cache plus cache summary. |
| `rc_build_target_microgem()` | Build one structural target-local micro-GEM. | Micro-GEM plus closure and medium diagnostics. |
| `rc_run_microcompass()` | Run cached structural micro-GEM strict LP with unit-specific penalties. | Score, penalty, vmax, feasibility, LP diagnostics. |
| `rc_test_microcompass_differential()` | Test sample-level differences with `lm`, `wilcox`, or continuous `limma::lmFit()`/`eBayes()`. | Differential result table. |
| `rc_export_microcompass()` | Write standardized outputs. | RDS/TSV files under numbered folders. |

## Required input tables

### Medium scenario table

Usually create this with `rc_make_medium_scenarios()`. A custom table must contain at least:

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

Common roles: `internal`, `boundary_like`, `exchange`, `transport`, `demand`, `sink`, `biomass`, `maintenance`, `cofactor_recycle`, `artificial_support`, `blocked`, `unknown`.

## Outputs to inspect

After `rc_export_microcompass(res, outdir)`, inspect:

```text
02_medium/medium_scenarios.tsv.gz
02_medium/medium_sensitivity_summary.tsv.gz
03_microgem/closure_diagnostics.tsv.gz
03_microgem/microgem_cache_summary.tsv.gz
04_microcompass/strict_score_matrix.rds
04_microcompass/strict_penalty_matrix.rds
04_microcompass/lp_diagnostics.tsv.gz
session_info.txt
```

## Interpretation

Use concise conservative language:

- `strict_feasible + low penalty + stable medium sensitivity`: stronger reaction-potential evidence.
- Strict infeasible LP results receive score 0 under the specified structural micro-GEM and medium scenario.
- ATAC confidence reflects chromatin-expression concordance, not regulatory causality.
- DNA methylation is not modeled.

## Legacy aliases

These aliases remain for compatibility: `rc_make_supercell2_metacells()`, `rc_import_supercell2_metacells()`, `rc_run_layer1_from_metacells()`, `rc_run_layer2_compass_lp()`, and `rc_layer2_apply_bounds()`.
