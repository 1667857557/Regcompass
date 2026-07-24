# Manual core reaction union-GEM scoring

`rc_regcompass_step_target_union()` is an optional step after Layer 1. It is intended for hypothesis-driven validation in which one or more reactions, or genes directly represented in GEM GPR rules, are selected as scoring targets.

The step performs the following operations:

1. Resolve the requested reaction IDs and gene-derived reactions.
2. Treat only those reactions as `score_target = TRUE`.
3. Expand the structural model through the core reactions' subsystems and shared KEGG, Reactome, and master-Rhea reaction identifiers.
4. Mark all annotation-expanded reactions as `model_only = TRUE`.
5. Add any further reactions required by the existing add-only FASTCORE feasibility completion.
6. Run directional microCOMPASS LPs on the resulting shared union GEM while returning scores only for the selected core reactions.

Model-only reactions are not emitted as reaction targets. When they have GPR evidence, their penalties remain part of the network-level LP objective; this prevents pathway context from becoming cost-free while keeping the output restricted to the requested cores.

## Select reactions directly

```r
selected <- rc_regcompass_step_target_union(
  layer1 = step4,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_target_union",
  core_reaction_ids = c("MAR04324", "MAR03964"),
  expansion_mode = "ordered_once",
  layer2_args = list(
    omega = 0.95,
    solver = "highs",
    time_limit = 60
  )
)
```

## Select reactions from genes

The default `gene_match = "complete_gpr"` only selects a reaction when the supplied genes cover at least one complete GPR AND group. This avoids treating one subunit of an obligate enzyme complex as a complete reaction core.

```r
selected <- rc_regcompass_step_target_union(
  layer1 = step4,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_glutathione_union",
  core_genes = c("GCLC", "GCLM", "GSS", "GSR", "G6PD", "PGD"),
  gene_match = "complete_gpr"
)
```

Use `gene_match = "any_direct"` only when intentionally selecting every reaction that directly contains at least one requested gene, including incomplete enzyme complexes.

## Inspect the union

```r
selected$global_core_reactions
selected$global_reaction_membership[
  , c(
    "reaction_id",
    "is_core",
    "score_target",
    "model_only",
    "inclusion_stage",
    "source_annotation"
  )
]
selected$summary
selected$microcompass$penalty
```

The step writes:

- `selected_core_reactions.tsv.gz`
- `union_reaction_membership.tsv.gz`
- `union_summary.tsv.gz`
- `scores/`
- `step_target_union.rds`

For several selected cores, all annotation-expanded memberships are deduplicated into one shared union GEM before scoring. `expansion_mode = "fixed_point"` permits transitive annotation expansion; `"ordered_once"` is the conservative default.
