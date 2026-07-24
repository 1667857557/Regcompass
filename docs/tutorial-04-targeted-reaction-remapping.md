# Tutorial Level 4: remap selected genes or reactions and run second-pass scoring

Use this tutorial after a completed **stepwise** `meta_module_gem` analysis when the original Layer 2 run has already scored all global core reactions, but a specific metabolic gene set or selected core reaction set should be used as anchors to map and score directly database-linked non-core reactions.

This is a targeted restart. It does **not** rerun GRN inference, metacell construction, meta-module construction, or Layer 1. It reuses the exact cached global union GEM and its medium-specific bounds.

## Required objects

```r
library(RegCompassR)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

The supplied objects must come from the same run. `step5` must have been generated with `model_mode = "meta_module_gem"`, and its model-cache files must still exist.

```r
stopifnot(
  inherits(step3, "regcompass_meta_module_step"),
  inherits(step4, "regcompass_layer1_step"),
  inherits(step5, "regcompass_layer2_step"),
  identical(step5$model_mode, "meta_module_gem"),
  identical(step3$gem_fingerprint, step4$gem_fingerprint),
  identical(step4$gem_fingerprint, step5$gem_fingerprint),
  all(file.exists(step5$model_cache_summary$file))
)
```

## Option A: select anchors by metabolic genes

`core_genes` is resolved through the GEM Boolean GPR table and then restricted to reactions that were global cores in the original Layer 2 analysis.

```r
targeted_by_gene <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/05b_targeted_GSH_PPP",
  core_genes = c(
    "SLC7A11",
    "GCLC",
    "GCLM",
    "GSS",
    "GSR",
    "G6PD",
    "PGD"
  ),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  ),
  parallel = TRUE,
  progress = TRUE
)
```

`gene_match = "complete_gpr"` requires all genes in at least one GPR AND-group to be present in the supplied gene set. Use `gene_match = "any_direct"` only when any direct GPR participation is intended; it is broader and can select multisubunit reactions for which the supplied genes alone are not sufficient.

Inspect how genes resolved to previous core anchors:

```r
targeted_by_gene$selected_core_reactions[
  ,
  c("gene", "reaction_id", "selection_source")
]
```

## Option B: select anchors by reaction IDs

Selected reaction IDs must already have been global core targets in the original Layer 2 run.

```r
available_core_ids <- unique(
  as.character(step3$global_modules$global_core_reactions$reaction_id)
)

head(available_core_ids)

selected_reaction_ids <- c("MAR04324", "MAR06231")
stopifnot(all(selected_reaction_ids %in% available_core_ids))

targeted_by_reaction <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/05b_targeted_reactions",
  core_reaction_ids = selected_reaction_ids,
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  ),
  parallel = TRUE,
  progress = TRUE
)
```

Genes and reaction IDs may be supplied together. The resolved anchor set is their union.

## What is remapped

For each selected previous core reaction, RegCompass directly identifies reactions sharing at least one:

- KEGG reaction identifier;
- Reactome reaction identifier;
- master-Rhea identifier.

Only mapped reactions that were **not** global core targets in the original Layer 2 run are scored in the second pass. The following are deliberately excluded:

- same-subsystem expansion by itself;
- recursive or transitive propagation;
- metabolite-neighbour expansion;
- FASTCORE-only support reactions without a direct database cross-reference;
- reactions already scored as original global cores.

Inspect the mapping relations and unique second-pass targets:

```r
targeted_by_gene$expanded_reaction_catalog[
  ,
  c(
    "anchor_core_reaction_id",
    "reaction_id",
    "expansion_type",
    "source_annotation",
    "previous_union_is_core",
    "score_target"
  )
]

targeted_by_gene$expanded_scoring_targets
```

Expected `expansion_type` values are:

```text
shared_kegg_reaction
shared_reactome_reaction
shared_master_rhea_reaction
```

## Confirm exact structural-model reuse

The second pass must reuse the exact cached union model rather than rebuild a new reaction network.

```r
stopifnot(
  inherits(targeted_by_gene, "regcompass_target_union_step"),
  all(targeted_by_gene$expanded_scoring_targets$score_target),
  targeted_by_gene$microcompass$params$structural_model_reused_exactly,
  identical(
    targeted_by_gene$microcompass$params$target_scope,
    "direct_kegg_reactome_master_rhea_noncore_only"
  ),
  all(
    targeted_by_gene$microcompass$model_cache_summary$
      reused_without_rebuilding
  )
)
```

## Interpret and export the second-pass scores

```r
dim(targeted_by_gene$microcompass$penalty)
dim(targeted_by_gene$microcompass$score)

targeted_by_gene$microcompass$lp_diagnostics
```

`penalty` remains the primary output: lower values indicate that the required target flux is more compatible with the integrated evidence. `score` is a within-target relative penalty rank, not a probability and not a measured flux.

The output directory contains:

```text
selected_previous_core_reactions.tsv.gz
expanded_reaction_catalog.tsv.gz
expanded_scoring_targets.tsv.gz
reused_global_union_membership.tsv.gz
target_union_summary.tsv.gz
scores/
step_target_union.rds
step_timing.tsv
```

Changing `core_genes` or `core_reaction_ids` only requires rerunning this tutorial. Changing GRNs, `gamma = 30`, meta-module construction, the GEM, medium bounds, or the original Layer 2 target model requires restarting from the earliest affected upstream stage.