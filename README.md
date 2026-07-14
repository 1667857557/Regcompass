# RegCompassR

RegCompassR 1.3 integrates sample-aware RNA+ATAC metacells, GPR-aware reaction evidence, sample-specific Pando regulatory networks, and steady-state reaction-potential analysis.

## Implemented workflow

```text
Seurat RNA+ATAC object
→ condition × sample × cell-type strata
→ SuperCell2 RNA+ATAC metacells
→ metacell fragment aggregation and LinkPeaks
→ Layer 1 RNA-GPR capacity + ATAC-supported confidence
→ sample-specific Pando metabolic GRNs
→ GRN genes mapped to Human-GEM core reactions
→ subsystem + KEGG/Reactome + master-Rhea biological expansion
→ full-GEM or FASTCORE-completed meta-module-GEM
→ directional microCOMPASS scoring
```

Layer 2 supports exactly two structural modes:

| Mode | Structural model |
|---|---|
| `full_gem` | The complete medium-constrained Human-GEM, analogous to the structural model used by COMPASS. |
| `meta_module_gem` | The complete GRN/subsystem/pathway biological reaction set plus compact support reactions selected from Human-GEM by add-only FASTCORE. |

Target-k-hop, module-meso-GEM, and automatic fallback modes were removed in version 1.3.

## Installation

```r
remotes::install_github("1667857557/SuperCell_Seurat_V4@supercell-2.0")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

The installed Pando package must retain repository metadata pointing to `1667857557/Pando_regcompass`.

## Prepare Human-GEM

```r
library(RegCompassR)

gem <- rc_prepare_human2_gem_v12(version = "2.0.0")
```

Human-GEM keeps one signed variable per reaction:

```text
lb < 0 < ub     reversible
0 <= lb < ub    forward-only
lb < ub <= 0    reverse-only
lb = ub = 0     blocked
```

RegCompass derives directions from numerical bounds. A metadata `reversible` flag does not override the active bounds.

## Medium constraints

For biological interpretation, provide a curated medium table rather than treating every exchange reaction as available.

```r
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

`rc_apply_medium_constraints()` first closes uptake through all annotated exchange reactions and then opens only the listed available exchanges. Built-in scenarios are broad sensitivity scenarios and should not be treated as tissue-specific media without additional curation.

## Integrated meta-module-GEM analysis

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
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000
    )
  )
)
```

The result is saved as:

```text
RegCompassR_v1.3_meta_module/regcompass_v1.3_result.rds
```

## Integrated full-GEM analysis

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

The integrated full-GEM workflow still uses the Pando-derived core-reaction list as its default scoring target set. Supply `target_reactions` inside `layer2_args` to score another explicit set.

## Biological meta-module definition

For each sample-specific connected metabolic GRN component, RegCompass performs the following ordered expansion:

1. Map GRN metabolic genes to Human-GEM GPR reactions. These reactions are marked `is_core = TRUE`.
2. Include all reactions assigned to each core reaction's subsystem.
3. Include reactions in subsystems linked by shared KEGG or Reactome reaction identifiers.
4. Include reactions in subsystems linked by shared master-Rhea identifiers.

```r
core_reactions <- result$grn_meta_modules$core_gene_reaction
reaction_membership <- result$grn_meta_modules$reaction_membership
```

The expanded set is a biological envelope. Every member is retained in the completed local GEM, but only direct GRN/GPR core reactions are mandatory directional FASTCORE targets. This prevents one blocked peripheral pathway annotation from invalidating the entire biological module.

## FASTCORE completion

A single completed model can be built directly:

```r
module_gem <- rc_build_meta_module_gem(
  gem = gem,
  reaction_membership = reaction_membership,
  core_reactions = core_reactions,
  sample_id = "sample_1",
  module_id = "sample_1::GRN0001",
  medium_table = medium,
  target_direction = "both",
  solver = "highs",
  fastcore_epsilon = 1e-4,
  strict = TRUE
)
```

The implementation performs:

1. Medium application to the complete Human-GEM.
2. Disabling of demand, sink, and artificial-support reactions for structural completion.
3. Parent-model feasibility validation.
4. FASTCC-style directional consistency screening.
5. Parent and biological-envelope directional checks for each core reaction.
6. FASTCORE LP-7 selection of simultaneously active unresolved core tasks.
7. FASTCORE LP-10 minimization of absolute flux outside the biological and already selected support sets.
8. Union of the biological envelope and selected support reactions.
9. Final validation of every parent-feasible target direction.

FASTCORE is add-only in RegCompass: it does not prune the biological reaction envelope.

### Reaction provenance

```r
module_gem$reaction_meta[, c(
  "reaction_id",
  "biological_meta_module_member",
  "fastcore_support",
  "support_only"
)]
```

### Structural diagnostics

```r
module_gem$closure_diagnostics
module_gem$completion_iterations
module_gem$build_params
module_gem$target_status
```

A parent-blocked reaction is reported as `parent_blocked`. RegCompass does not create an artificial exchange, demand, sink, or gap-filling reaction to make it feasible.

## Direct Layer 2 analysis

### Full GEM

```r
full_result <- rc_run_microcompass(
  layer1 = layer1,
  gem = gem,
  target_reactions = target_reactions,
  medium_scenarios = medium,
  mode = "full_gem",
  unit = "sample_celltype",
  target_direction = "both",
  solver = "highs"
)
```

### Meta-module GEM

```r
module_result <- rc_run_microcompass(
  layer1 = layer1,
  gem = gem,
  target_reactions = unique(core_reactions$reaction_id),
  medium_scenarios = medium,
  mode = "meta_module_gem",
  reaction_membership = reaction_membership,
  core_reactions = core_reactions,
  unit = "sample_celltype",
  target_direction = "both",
  solver = "highs",
  model_params = list(
    fastcore_epsilon = 1e-4,
    strict = TRUE
  )
)
```

A completed meta-module is cached once per:

```text
sample_id × module_id × medium_scenario
```

All core reactions from that module reuse the same structural model. The model is evaluated only against Layer 2 units carrying the matching `sample_id`.

## Directional microCOMPASS mathematics

For direction sign `d` (`+1` forward and `-1` reverse), step 1 computes:

```text
maximize d × v_target
subject to S v = 0
           lb <= v <= ub
```

Step 2 minimizes evidence-weighted absolute flux:

```text
minimize sum_i penalty_i × a_i
subject to S v = 0
           lb <= v <= ub
           -a_i <= v_i <= a_i
           a_i >= 0
           d × v_target >= omega × vmax
```

The signed-flux formulation preserves forced non-zero positive or negative reaction bounds. In meta-module mode, `vmax` is always calculated in the completed local GEM itself; there is no alternative reference-mode parameter.

## Main output fields

```r
result$microcompass$score
result$microcompass$penalty
result$microcompass$vmax
result$microcompass$feasible
result$microcompass$evaluated
result$microcompass$model_cache_summary
result$microcompass$model_diagnostics
result$microcompass$lp_diagnostics
```

In meta-module mode, score matrix row IDs include sample and module identity:

```text
sample=<sample>::module=<module>::reaction=<reaction>::direction=<direction>::medium=<medium>
```

## Interpretation limits

- FASTCORE is an LP-based compact reconstruction algorithm, not an exact minimum-cardinality MILP.
- Steady-state feasibility does not establish thermodynamic or kinetic feasibility.
- Reaction capacity and microCOMPASS scores are multiome-supported reaction potentials, not measured fluxes.
- The completed model guarantees requested parent-feasible core directions; it does not require every peripheral biological-envelope reaction to carry flux.
- Results remain conditional on Human-GEM stoichiometry, bounds, reaction-role annotations, GPR mapping, and the selected medium.

See `docs/meta_module_v13_design.md` for the mathematical and engineering specification.
