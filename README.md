# RegCompassR workflow

RegCompassR 1.2 runs a strict multiome workflow with an optional sample-specific Pando GRN meta-module layer. The implemented biological-envelope stages and the intended local-GEM support-completion stages are:

```text
Seurat RNA+ATAC object
→ pre-filter condition × sample × cell type strata
→ build SuperCell2 RNA+ATAC metacells inside each retained stratum
→ post-filter strata by actual metacell count
→ aggregate ATAC fragments separately for each retained stratum
→ recompute stratum-specific metacell LinkPeaks
→ Layer 1 RNA+ATAC-GPR capacity/confidence
→ sample-specific Pando GRN
→ GRN core reactions
→ subsystem + KEGG/Reactome + master-Rhea expansion
→ biological reaction envelope B
→ medium-specific full GEM
→ FASTCC consistent parent
→ core direction precheck
→ add-only FASTCORE support completion
→ capacity-retention validation
→ completed local GEM = B* ∪ support
→ sample/module-aware microCOMPASS
```

The Pando/meta-module layer defines a biological reaction envelope **B** before solver support is considered. In the current code, **B** is implemented by Layer 1 RNA+ATAC-GPR evidence, sample-specific Pando GRNs, GRN core reactions, and subsystem + KEGG/Reactome + master-Rhea expansion. The intended complete local-GEM path after **B** is medium-specific full GEM → FASTCC consistent parent → core-direction precheck → add-only FASTCORE support completion → capacity-retention validation → completed local GEM (**B*** ∪ support) → sample/module-aware microCOMPASS. The current `rc_build_meta_module_gem()` implementation has not yet added FASTCC/FASTCORE/capacity-retention steps; it currently performs deterministic Human-GEM support inclusion and keeps those support reactions separate from biological members.

## Install the pinned dependencies

```r
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

The Pando package name is `Pando`. Runtime validation checks that the installed package metadata points to `1667857557/Pando_regcompass`; it records but does not restrict the package version or commit SHA.

## RegCompassR 1.2 example

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

gem <- rc_prepare_human2_gem_v12(version = "2.0.0")

result <- rc_run_regcompass_v12(
  object = object,
  gem = gem,
  outdir = "RegCompassR_v1.2_run",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  rna_assay = "RNA",
  atac_assay = "ATAC",
  metacell_args = list(
    gamma = 150,
    adaptive_gamma = TRUE,
    min_cells_pre_metacell = 100,
    min_metacell_size = 10,
    min_metacells_post_metacell = 10,
    future_plan = "sequential"
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
    min_shared_tfs = 1,
    expansion_mode = "ordered_once"
  )
)

meta_modules <- result$grn_meta_modules$reaction_membership
module_summary <- result$grn_meta_modules$meta_module_summary
```

The metabolic Pando target genes are the case-insensitive intersection of original single-cell RNA features, retained metacell RNA features, and Human-GEM GPR genes. Pando is run separately for each `sample_id`.


## Function-by-function workflow guide

The v1.2 API is intentionally split into auditable functions. Use them in this order when you need to inspect or rerun individual stages instead of calling `rc_run_regcompass_v12()` as a wrapper:

| Step | Function | Actual implementation role | Main outputs / side effects |
|---|---|---|---|
| 1 | `rc_prepare_human2_gem_v12()` | Downloads or reads a pinned Human-GEM release, converts `Human-GEM.yml`, enriches reaction metadata with subsystem, KEGG, Reactome, Rhea and master-Rhea IDs, annotates reaction roles, and stores GPR/metabolic gene tables. | A validated GEM list with `gpr_table`, `metabolic_genes`, `reaction_meta`, `reaction_roles`, and `model_info$annotation_schema = "regcompass_humangem_v12"`. |
| 2 | `rc_run_regcompass_multiome_metacell()` | Runs the existing strict Layer 1 workflow: pre-filter strict strata, build SuperCell2 RNA+ATAC metacells, post-filter by actual metacell count, aggregate fragments, recompute stratum-specific LinkPeaks, then compute GPR capacity and ATAC confidence. | Layer 1 result plus saved `00_stratum_qc/`, `01_metacells/`, `02_metacell_fragments/`, and `03_linkpeaks/` outputs. |
| 3 | `rc_load_metacell_object_from_run()` | Reloads saved `metacell_object.rds` files from `01_metacells/`, optionally subsets to retained metacell IDs, clears Signac fragment pointers, and merges sample/stratum objects. | One retained metacell Seurat object for the Pando stage. |
| 4 | `rc_run_pando_meta_modules()` | Validates the pinned Pando fork, derives metabolic target genes from original RNA genes ∩ retained metacell RNA genes ∩ Human-GEM metabolic genes, runs Pando once per sample, records per-sample status/errors, projects significant Pando edges to metabolic gene modules, maps modules to reactions, expands reactions, and writes audit tables. | `04_pando_meta_modules/` files and a `regcompass_pando_meta_module_v1.2` list. |
| 4a | `rc_extract_pando_tf_peak_gene()` | Converts `coef(grn_object)` to a table, merges `Pando::gof()` by target when available, uppercases TF/target names, and filters by `padj`, `estimate`, and optional `rsq`. | Complete and significant TF–peak–gene tables. |
| 4b | `rc_project_metabolic_grn()` | Keeps significant metabolic targets, creates shared-TF gene pairs, optionally adds direct metabolic-TF edges, applies shared-TF/Jaccard/top-k filters, and labels connected components as sample-specific `GRN####` modules. | `metabolic_gene_nodes` and `metabolic_gene_edges`. |
| 4c | `rc_map_meta_module_core_reactions()` | Uppercase-joins GRN module genes to Human-GEM GPR genes and marks every matching reaction as a core GRN reaction. | `core_gene_reaction` with `inclusion_stage = "core_grn_gene"`. |
| 4d | `rc_expand_meta_module_reactions()` | For each `sample_id × module_id`, filters to valid GEM reactions, adds core-subsystem reactions, then KEGG/Reactome-linked subsystem reactions, then master-Rhea-linked subsystem reactions; `fixed_point` repeats until stable or `max_iterations`. | `reaction_membership`, `meta_module_summary`, and normalized cross-reference maps. |
| 5 | Intended FASTCC/FASTCORE local-GEM completion | After biological envelope **B**, the intended flow is medium-specific full GEM, FASTCC consistent parent, core-direction precheck, add-only FASTCORE support completion, and capacity-retention validation. | Completed local GEM = **B*** ∪ support, with support reactions kept outside biological membership. |
| 5 current | `rc_build_meta_module_gem()` | Current implementation uses a selected biological meta-module as `grn_meta_module`, delegates deterministic local closure/support expansion to `rc_build_module_meso_gem()`, and labels biological versus support-only reactions. | A local module GEM with `biological_meta_module_member` and `support_only` flags; FASTCC/FASTCORE/capacity-retention are not yet implemented here. |
| 6 | `rc_run_regcompass_v12()` | Convenience wrapper that runs Step 2, reloads retained metacells, runs Step 4, attaches `grn_meta_modules`, sets `schema_version`, and saves `regcompass_v1.2_result.rds`. | Integrated Layer 1 + GRN meta-module result. |

A staged tutorial mirroring these functions is maintained in `docs/meta_module_v12_design.md`.

## Meta-module definition

For each sample-specific connected metabolic GRN component:

1. Map all GRN metabolic genes to Human-GEM GPR reactions; these are core reactions.
2. Include every reaction from every subsystem assigned to a core reaction.
3. Collect KEGG/Reactome reaction identifiers from the current subsystem collection, identify all subsystems containing a matching identifier, and include all reactions in those subsystems.
4. Collect master-Rhea identifiers from all currently included reactions, identify all subsystems containing a reaction with a matching master-Rhea identifier, and include all reactions in those subsystems.

The default `expansion_mode = "ordered_once"` applies the rules once. `"fixed_point"` repeatedly applies them until stable and is intended for sensitivity analysis because it can produce much larger modules.

`UNASSIGNED`, `NA`, and `NONE` are not treated as valid shared subsystem labels.

## Building a local meta-module GEM

```r
module_gem <- rc_build_meta_module_gem(
  gem = gem,
  reaction_membership = meta_modules,
  sample_id = module_summary$sample_id[[1]],
  module_id = module_summary$module_id[[1]],
  medium_table = medium,
  include_one_hop = FALSE,
  include_transport = TRUE,
  include_exchange = TRUE,
  include_protected = TRUE
)
```

The exact expanded biological meta-module reaction set (**B**) is retained as `biological_meta_module_member = TRUE`. Reactions added only for transport, exchange, local closure, or solver feasibility are labeled `support_only = TRUE`; they must not redefine **B**. In the intended FASTCORE-based completion, FASTCORE is add-only and only contributes the minimal support set needed for core/biological reactions to carry steady-state flux after medium-specific full-GEM consistency and direction prechecks.

See `docs/meta_module_v12_design.md` for the algorithm, thresholds, output schema, and validation requirements.

## Existing Layer 1 and Layer 2 workflow

The existing formal entry point remains available:

```r
library(RegCompassR)

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
  BPPARAM_metacell = FALSE,
  BPPARAM_linkpeaks = bp,
  BPPARAM_layer1 = bp
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

## Strict metacell strata and LinkPeaks gates

The formal metacell workflow uses one fixed stratum definition for metacell construction, post-metacell filtering, fragment aggregation, and LinkPeaks:

```text
condition × sample × cell type
```

The workflow applies two auditable hard filters:

1. Before SuperCell, each stratum must contain at least `min_cells_pre_metacell` original cells.
2. After metacell construction, the same stratum must contain at least `min_metacells_post_metacell` generated metacells.
3. ATAC fragments are aggregated independently for each retained strict stratum.
4. LinkPeaks is recomputed independently in each retained strict stratum.
5. Excluded cells and metacells remain only in QC reports under `00_stratum_qc/`.

`state_col` can remain in metadata for later summaries but does not alter the formal stratum.

## Layer 2 GEM cache strategies

`rc_run_microcompass()` defaults to a full-GEM analysis. Alternative `microgem_params$strategy` values are:

| Strategy | Behavior |
|---|---|
| `"full_gem"` | Cache the complete GEM per medium scenario. |
| `"module_meso_gem"` | Cache one module-level meso-GEM per module and medium scenario. |
| `"target_khop"` | Build target-local k-hop micro-GEMs. |
| `"auto"` | Try target-local micro-GEMs and fall back to module meso-GEMs when closure fails. |

## RegCompassR 1.2 outputs

```text
04_pando_meta_modules/
├── pando_sample_status.tsv.gz
├── pando_tf_peak_gene_all.tsv.gz
├── pando_tf_peak_gene_significant.tsv.gz
├── metabolic_gene_nodes.tsv.gz
├── metabolic_gene_edges.tsv.gz
├── core_gene_reaction.tsv.gz
├── meta_module_reactions.tsv.gz
├── meta_module_summary.tsv.gz
├── pando_meta_modules.rds
├── sample_metacell_objects/<sample>.rds
└── pando_objects/<sample>.rds
```

The primary audit table is `meta_module_reactions.tsv.gz`, which records sample, GRN module, reaction, core status, inclusion stage, source annotations, and expansion mode.
