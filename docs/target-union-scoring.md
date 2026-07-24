# Direct database-linked non-core scoring in the existing union GEM

`rc_regcompass_step_target_union()` is an optional analysis after a completed stepwise Layer 2 run with `model_mode = "meta_module_gem"`.

## Exact target definition

A selected core reaction is a **mapping anchor**, not a second-pass LP target. Candidate reactions are obtained directly from identifiers attached to that core:

- shared KEGG reaction ID;
- shared Reactome reaction ID;
- shared master-Rhea ID.

No same-subsystem expansion is used. No candidate reaction is used as a new anchor, so the mapping is not transitive or recursive.

The second LP pass scores only directly linked reactions that were **not global core targets in the original Layer 2 run**. It excludes:

- the selected core itself;
- any other global core already scored in Layer 2;
- reactions related only by subsystem;
- reactions reachable only through an intermediate mapped reaction;
- FASTCORE-only support reactions without a direct database link to the selected core;
- generic union-GEM members and metabolite-neighbour reactions.

## Processing sequence

1. Resolve previous core anchors from reaction IDs or GPR genes.
2. Read each anchor's KEGG, Reactome, and master-Rhea identifiers.
3. Find reactions that directly share those identifiers with that anchor.
4. Preserve one relation row per `anchor × reaction × database mapping`.
5. Remove reactions already scored as global cores.
6. Require every remaining target to exist in the original global union GEM.
7. Reuse the exact cached stoichiometry, bounds, and medium.
8. Run directional Vmax/minimum-penalty LP only for the remaining non-core targets.

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

## Inspect direct mappings and LP targets

```r
expanded$selected_core_reactions

expanded$expanded_reaction_catalog[, c(
  "anchor_core_reaction_id",
  "reaction_id",
  "expansion_type",
  "source_annotation",
  "previous_union_is_core",
  "score_target",
  "lp_exclusion_reason"
)]

expanded$expanded_scoring_targets[, c(
  "reaction_id",
  "anchor_core_reaction_ids",
  "expansion_types",
  "source_annotations"
)]
```

Valid `expansion_type` values are:

```r
c(
  "shared_kegg_reaction",
  "shared_reactome_reaction",
  "shared_master_rhea_reaction"
)
```

A directly linked reaction that was already a global core remains in `expanded_reaction_catalog` with `score_target = FALSE` and `lp_exclusion_reason = "already_scored_in_original_layer2"`.

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
