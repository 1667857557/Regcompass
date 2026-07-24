# Re-score annotation-related reactions in the existing union GEM

`rc_regcompass_step_target_union()` is an optional **post-Layer-2** step. The normal RegCompass workflow must first finish the original LP analysis of all GRN-derived core reactions and persist its global union-GEM cache.

The second pass then performs the following operations:

1. Select one or more reactions that were core LP targets in the previous analysis, either by reaction ID or by their directly associated GPR genes.
2. Expand those selected cores through:
   - the same subsystem;
   - shared KEGG reaction IDs;
   - shared Reactome reaction IDs;
   - shared master-Rhea IDs.
3. Verify that every expanded reaction is already present in the previously constructed global union GEM.
4. Load the exact union-GEM model file recorded by the first Layer-2 run for each medium scenario.
5. Run a second directional LP pass in which **the selected cores and every annotation-related reaction are all scoring targets**.

The second pass does not rebuild the union GEM and does not classify expanded reactions as model-only support. Its purpose is to obtain LP scores for the wider pathway/function context surrounding selected core reactions while keeping the same stoichiometric network and medium-specific bounds used by the original core analysis.

## Required preceding steps

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules"
)

step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1"
)

step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_core_layer2",
  model_mode = "meta_module_gem"
)
```

`step5$model_cache_summary$file` must remain available because these files contain the previously constructed medium-specific global union GEMs.

## Select previously scored core reactions

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
    omega = 0.95,
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  )
)
```

A reaction supplied through `core_reaction_ids` must already occur in `step3$global_modules$global_core_reactions`; the function rejects arbitrary non-core reactions.

## Select previous core reactions from genes

The default `gene_match = "complete_gpr"` only resolves a previous core reaction when the supplied genes cover at least one complete GPR AND group. This avoids treating one subunit of an obligate enzyme complex as a complete core reaction.

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

Use `gene_match = "any_direct"` only when intentionally allowing a gene to resolve a previous core reaction despite incomplete coverage of an enzyme complex.

## Inspect the second-pass targets and scores

```r
expanded$selected_core_reactions
expanded$expanded_scoring_targets[
  , c(
    "reaction_id",
    "selected_core_anchor",
    "score_target",
    "target_role",
    "inclusion_stage",
    "source_annotation",
    "previous_union_inclusion_stage"
  )
]
expanded$summary
expanded$microcompass$penalty
expanded$microcompass$feasible
```

Every row in `expanded_scoring_targets` has `score_target = TRUE`. The LP output therefore contains the original selected cores plus all same-subsystem, KEGG/Reactome-linked, and master-Rhea-linked reactions that were present in the original union GEM.

## Outputs

The step writes:

- `selected_previous_core_reactions.tsv.gz`
- `expanded_scoring_targets.tsv.gz`
- `reused_global_union_membership.tsv.gz`
- `target_union_summary.tsv.gz`
- `scores/`
- `step_target_union.rds`

`expansion_mode = "ordered_once"` is the conservative default and matches the canonical annotation-expansion order. `"fixed_point"` is allowed only when all transitively expanded reactions already exist in the previously built global union GEM; otherwise the function stops rather than silently rebuilding or changing the model.
