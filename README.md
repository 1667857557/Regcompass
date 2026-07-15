# RegCompassR

RegCompassR 1.3 runs a sample-aware RNA+ATAC workflow that converts metacells, GPR evidence, sample-specific Pando regulatory modules, Human-GEM structure, and medium constraints into directional microCOMPASS reaction-potential scores.

## Main workflow

```text
Seurat RNA+ATAC object
→ condition × sample × cell-type strata
→ SuperCell2 RNA+ATAC metacells and fragment aggregation
→ LinkPeaks relinking within retained strata
→ Layer 1 RNA-GPR capacity and ATAC-supported confidence
→ sample-specific Pando metabolic GRNs
→ GRN genes mapped to Human-GEM core reactions
→ subsystem + KEGG/Reactome + master-Rhea meta-module expansion
→ full-GEM or FASTCORE-completed meta-module-GEM
→ directional microCOMPASS scoring
```

RegCompass 1.3 intentionally exposes two structural model modes:

| Mode | Structural model |
|---|---|
| `meta_module_gem` | The GRN/subsystem/pathway biological reaction set completed with compact Human-GEM support reactions by add-only FASTCORE. This is the recommended default. |
| `full_gem` | The complete medium-constrained Human-GEM, analogous to the structural model used by COMPASS. |

Older target-k-hop, module-meso-GEM, automatic-fallback, and alternate reference-mode APIs are not part of the current tutorial.

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

`rc_apply_medium_constraints()` first closes uptake through all annotated exchange reactions and then opens only the listed available exchanges. Built-in medium scenarios are broad sensitivity scenarios; use curated tables for biological interpretation.

Human-GEM reaction directions are derived from active numerical bounds:

```text
lb < 0 < ub     reversible
0 <= lb < ub    forward-only
lb < ub <= 0    reverse-only
lb = ub = 0     blocked
```

## Recommended integrated analysis

Use `rc_run_regcompass()` for the supported end-to-end workflow. The default tutorial path uses `model_mode = "meta_module_gem"`.

Parallelization is split into two stages: metacell construction and Pando/meta-module inference run by retained `condition × sample × cell-type` strata, while the downstream COMPASS-like Layer 2 task grid uses its own parallel pass after the upstream worker pool has been released. By default, sample-derived meta-modules are merged into one completed meta-module GEM so every metacell is scored on the same structure; metacell-specific RNA+ATAC evidence still enters through each unit's own penalty vector.

```r
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

result <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "RegCompassR_v1.3_meta_module",
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
  ),
  layer2_args = list(
    unit = "sample_celltype",
    target_direction = "both",
    solver = "highs",
    parallel = TRUE,
    model_params = list(
      strict = TRUE,
      merge_sample_modules = TRUE,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000
    )
  )
)
```

The result is saved to:

```text
RegCompassR_v1.3_meta_module/regcompass_v1.3_result.rds
```

## Optional full-GEM analysis

Switch only `model_mode` and the output directory when you want to score targets inside the full Human-GEM structure.

```r
result_full <- rc_run_regcompass(
  object = object,
  gem = gem,
  outdir = "RegCompassR_v1.3_full_gem",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = fragment_files,
  model_mode = "full_gem",
  medium_scenarios = medium,
  layer2_args = list(
    unit = "sample_celltype",
    target_direction = "both",
    solver = "highs"
  )
)
```

The full-GEM workflow still uses Pando-derived core reactions as the default target set. Provide `target_reactions` inside `layer2_args` only when scoring an explicit alternative set.

## Key outputs

```r
result$grn_meta_modules$core_gene_reaction
result$grn_meta_modules$reaction_membership
result$microcompass$score
result$microcompass$penalty
result$microcompass$vmax
result$microcompass$feasible
result$microcompass$evaluated
result$microcompass$model_cache_summary
result$microcompass$model_diagnostics
result$microcompass$lp_diagnostics
```

In meta-module mode, score row IDs include sample and module identity:

```text
sample=<sample>::module=<module>::reaction=<reaction>::direction=<direction>::medium=<medium>
```

## Interpretation limits

- FASTCORE is an LP-based compact reconstruction algorithm, not an exact minimum-cardinality MILP.
- Steady-state feasibility does not establish thermodynamic or kinetic feasibility.
- Reaction capacity and microCOMPASS scores are multiome-supported reaction potentials, not measured fluxes.
- The completed model guarantees requested parent-feasible core directions; it does not require every peripheral biological-envelope reaction to carry flux.
- Results remain conditional on Human-GEM stoichiometry, bounds, reaction-role annotations, GPR mapping, and the selected medium.
- Metacells are technical aggregation units and are never treated as independent biological replicates.

See [`docs/meta_module_v13_design.md`](docs/meta_module_v13_design.md) for the engineering specification and [`docs/math_biology_audit_v13.md`](docs/math_biology_audit_v13.md) for the LP equations, biological assumptions, and validation notes.
