# Expanded non-core target scoring in the existing union GEM

`rc_regcompass_step_target_union()` is an optional analysis after a completed stepwise Layer 2 run with `model_mode = "meta_module_gem"`.

## Exact target definition

A selected core reaction is an **expansion anchor**, not a second-pass LP target. The function derives candidate reactions only through:

- the same subsystem as the selected core;
- shared KEGG reaction identifiers;
- shared Reactome reaction identifiers;
- shared master-Rhea identifiers.

The second LP pass scores only those annotation-expanded reactions that were **not global core targets in the original Layer 2 run**. It excludes:

- the selected core reaction itself;
- any other global core reaction already scored in Layer 2;
- local FASTCORE support reactions unless they independently qualify through one of the four annotation mappings above;
- generic union-GEM members;
- metabolite-neighbour or one-hop reactions.

## Processing sequence

1. Select previous core anchors by reaction ID or GPR gene.
2. Expand the anchors through subsystem, KEGG, Reactome, and master-Rhea mappings.
3. Retain the complete expansion catalog for provenance.
4. Remove reactions already scored as global cores.
5. Require every remaining target to exist in the original global union GEM.
6. Reuse the exact cached stoichiometry, bounds, and medium for each scenario.
7. Run the standard directional Vmax/minimum-penalty LP only for the remaining non-core expansion targets.

The function does not rerun FASTCORE or rebuild the union GEM.

## Required inputs

```r
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

The files listed in `step5$model_cache_summary$file` must remain available. Stage classes, workflow parameters, GEM fingerprints, metacell order, the original core set, and source-model hashes are checked before scoring.

## Select core anchors by reaction ID

```r
expanded <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_expanded_target_scoring",
  core_reaction_ids = c("MAR04324", "MAR03964"),
  expansion_mode = "ordered_once",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    omega = 0.95,
    time_limit = 60
  )
)
```

Reaction IDs must be members of `step3$global_modules$global_core_reactions` and must have been scored by `step5`.

## Select core anchors by gene

```r
expanded <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_glutathione_context",
  core_genes = c("GCLC", "GCLM", "GSS", "GSR", "G6PD", "PGD"),
  gene_match = "complete_gpr"
)
```

`complete_gpr` requires the supplied genes to cover one complete GPR AND group. Use `any_direct` only when partial enzyme-complex matching is intentional.

## Inspect the target split

```r
expanded$selected_core_reactions

expanded$expanded_reaction_catalog[, c(
  "reaction_id",
  "selected_core_anchor",
  "previous_union_is_core",
  "inclusion_stage",
  "source_annotation",
  "score_target",
  "lp_exclusion_reason"
)]

expanded$expanded_scoring_targets[, c(
  "reaction_id",
  "inclusion_stage",
  "source_annotation"
)]
```

Valid scoring-target stages are:

```r
c(
  "same_core_subsystem",
  "shared_kegg_or_reactome_reaction",
  "shared_master_rhea_reaction"
)
```

Core rows remain in `expanded_reaction_catalog` with `score_target = FALSE` and `lp_exclusion_reason = "already_scored_in_original_layer2"`.

```r
expanded$microcompass$penalty
expanded$microcompass$feasible
expanded$microcompass$model_cache_summary[, c(
  "medium_scenario",
  "file",
  "source_model_fingerprint",
  "source_model_md5",
  "reused_without_rebuilding"
)]
```

The raw minimum penalty is the primary result. Lower values indicate stronger compatibility with the integrated evidence under the fixed union-GEM constraints.

## Files

- `selected_previous_core_reactions.tsv.gz`
- `expanded_reaction_catalog.tsv.gz`
- `expanded_scoring_targets.tsv.gz`
- `reused_global_union_membership.tsv.gz`
- `target_union_summary.tsv.gz`
- `scores/`
- `step_target_union.rds`

`ordered_once` is the conservative expansion mode. `fixed_point` is allowed only when every transitively annotation-linked target already exists in the original union GEM.
