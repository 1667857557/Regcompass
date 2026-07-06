# RegCompassR microCOMPASS quick guide

RegCompassR analyzes target-local, medium-aware, multiome-supported reaction potential from annotated Seurat/Signac RNA+ATAC data. It does not infer whole-GEM fluxes, true uptake/secretion fluxes, or enzyme activity.

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
  outdir = "RegCompassR_run/01_metacells",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  rna_assay = "RNA",
  atac_assay = "ATAC",
  rna_reduction = "pca",
  atac_reduction = "lsi",
  gamma = 75,
  min_cells_per_stratum = 100,
  min_metacell_size = 20,
  filter_low_power_metacells = TRUE
)

gem <- rc_read_gem("HumanGEM_regcompass.rds")

layer1 <- rc_run_layer1_multiome(
  gpr_table = gem$gpr_table,
  rna_metacell_counts = mc$rna_counts,
  metacell_meta = mc$metacell_meta_used,
  atac_metacell_counts = mc$atac_counts,
  peak_gene_links = peak_gene_links,
  stratum_col = "cell_type"
)

gem <- rc_annotate_reaction_roles(gem, reaction_role_table = reaction_roles)
gem <- rc_apply_medium_constraints(gem, medium_table = medium)$gem

targets <- rc_select_target_reactions(
  layer1,
  method = "top_capacity",
  top_n = 100,
  min_C_rel = 0.15,
  min_confidence = 0.25
)

res <- rc_run_microcompass(
  layer1 = layer1,
  gem = gem,
  target_reactions = targets$reaction_id,
  medium_table = medium,
  unit = "sample_celltype",
  target_direction = "both",
  run_relaxed = TRUE,
  run_fva = TRUE,
  solver = "highs"
)

stat <- rc_test_microcompass_differential(
  res,
  formula = score ~ condition,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition"
)

rc_export_microcompass(res, "RegCompassR_run")
```

## Function checklist

| Step | Function | What it expects | What it returns |
|---|---|---|---|
| Validate input | `rc_validate_multiome_input()` | Seurat object with raw RNA counts, optional ATAC counts, and metadata columns | Invisible `TRUE`; stops on missing assays, metadata, or RNA/ATAC barcode mismatch |
| Build metacells | `rc_make_metacells()` | Annotated Seurat object and SuperCell reductions | Metacell counts, metadata, membership, diagnostics, `metacell_meta_all`, `metacell_meta_used` |
| Import metacells | `rc_import_metacells()` | Saved metacell directories | Same count/metadata structure as `rc_make_metacells()` |
| Run Layer 1 | `rc_run_layer1_multiome()` | GPR table, raw metacell RNA counts, metacell metadata, optional ATAC/link evidence | `C_rel`, confidence tables, diagnostics, sample-level Layer 1 matrices where available |
| Read GEM | `rc_read_gem()` | RDS file containing a GEM list | Validated GEM list |
| Validate GEM | `rc_validate_gem()` | GEM list with sparse/dense `S`, `lb`, `ub` | Sparse `S`, aligned bounds, IDs, and validation diagnostics |
| Annotate roles | `rc_annotate_reaction_roles()` | GEM and optional curated role table | GEM with `reaction_meta$role`, `role_source`, and `role_confidence` |
| Apply medium | `rc_apply_medium_constraints()` | GEM and medium table | `list(gem, medium_diagnostics)` |
| Select targets | `rc_select_target_reactions()` | Layer 1 result | Target reaction table; it does not expand networks |
| Build micro-GEM | `rc_build_target_microgem()` | GEM, one target reaction, optional medium table | Target-local GEM plus closure, medium, and gapfill diagnostics |
| Check closure | `rc_check_microgem_closure()` | Micro-GEM and target reaction | Strict target feasibility and boundary/dead-end counts |
| Run microCOMPASS | `rc_run_microcompass()` | Layer 1, GEM, target IDs | Score, penalty, vmax, feasibility, LP diagnostics, optional relaxed/FVA outputs |
| Relaxed LP | `rc_run_relaxed_balance_lp()` | Micro-GEM, penalties, target reaction | Slack feasibility and top slack diagnostics |
| Selected FVA | `rc_run_selected_fva()` | Micro-GEM and target reaction | Selected reaction min/max ranges and blocked/alternative-route flags |
| Differential test | `rc_test_microcompass_differential()` | `rc_run_microcompass()` result | Sample-level condition test table |
| Export | `rc_export_microcompass()` | microCOMPASS result and output directory | RDS and TSV diagnostics under standardized folders |

## Required table formats

### Medium table

```text
exchange_reaction_id  metabolite_id  condition  lb   ub    available
EX_glc_D_e            glc_D_e        all       -10  1000  TRUE
EX_gln_L_e            gln_L_e        all        -5  1000  TRUE
EX_lac_L_e            lac_L_e        all         0  1000  TRUE
```

### Reaction role table

```text
reaction_id  role       role_source
EX_glc_D_e   exchange   curated
R_HEX1       internal   curated
```

Recognized roles include `internal`, `exchange`, `transport`, `demand`, `sink`, `biomass`, `maintenance`, `cofactor_recycle`, `artificial_support`, `blocked`, and `unknown`.

## Result interpretation

Use conservative language:

- OK: target-local multiome-supported reaction potential differs by condition.
- OK: the target is feasible in the target micro-GEM under the specified medium.
- OK: high slack or wide FVA means the local-network interpretation is weak.
- Not OK: the result proves true flux, secretion, uptake, enzyme activity, or ATAC causality.

## Legacy names

The older `rc_make_supercell2_metacells()`, `rc_import_supercell2_metacells()`, `rc_run_layer1_from_metacells()`, `rc_run_layer2_compass_lp()`, and `rc_layer2_apply_bounds()` names remain for compatibility. Prefer the shorter microCOMPASS names shown above.
