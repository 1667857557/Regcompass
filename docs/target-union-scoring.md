# Expanded target scoring in the existing union GEM

`rc_regcompass_step_target_union()` is an optional post-Layer-2 analysis. It requires a completed stepwise run with `model_mode = "meta_module_gem"`.

## What it does

1. Select previous core targets by reaction ID or GPR gene.
2. Expand those anchors to reactions in the same subsystem or sharing KEGG, Reactome, or master-Rhea identifiers.
3. Require every expanded reaction to already exist in the original global union GEM.
4. Reuse the exact cached model file, stoichiometry, bounds, and medium for each scenario.
5. Run the standard directional Vmax/minimum-penalty LP for every expanded reaction.

The function does not rerun FASTCORE, rebuild the union, or classify expanded reactions as model-only support.

## Required inputs

```r
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

The files listed in `step5$model_cache_summary$file` must remain available. The function verifies stage classes, workflow parameters, GEM fingerprints, metacell order, the original core set, and source-model hashes before scoring.

## Select core reactions

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

Reaction IDs must be members of `step3$global_modules$global_core_reactions` and must have been used by `step5`.

## Select core reactions by gene

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

## Inspect results

```r
expanded$selected_core_reactions
expanded$expanded_scoring_targets[, c(
  "reaction_id",
  "selected_core_anchor",
  "target_role",
  "inclusion_stage",
  "source_annotation",
  "previous_union_inclusion_stage"
)]

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

Every row in `expanded_scoring_targets` is an LP target. The raw minimum penalty is the primary result; lower values indicate stronger compatibility with the integrated evidence under the fixed union-GEM constraints.

## Files

- `selected_previous_core_reactions.tsv.gz`
- `expanded_scoring_targets.tsv.gz`
- `reused_global_union_membership.tsv.gz`
- `target_union_summary.tsv.gz`
- `scores/`
- `step_target_union.rds`

`ordered_once` is the conservative expansion mode. `fixed_point` is allowed only when every transitively linked reaction is already present in the original union GEM.
