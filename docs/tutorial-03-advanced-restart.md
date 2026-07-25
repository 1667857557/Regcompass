# Tutorial Level 3: restart, sensitivity, and diagnostics

Use saved classed stage objects. RegCompass rejects cross-run object mixing when workflow settings, GEM fingerprints, metacell provenance, merged core sets, or ordered scoring units differ.

## Load a completed stepwise run

```r
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_steps/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

## Restart boundaries

### Rerun Stage 2 onward

Rerun metacells and every downstream stage after changing:

- `gamma`;
- RNA PCA/Harmony reduction or dimensions;
- ATAC LSI reduction or dimensions;
- cell membership or condition/cell-type metadata;
- RNA or ATAC count matrices.

### Rerun Stage 3 onward

Rerun meta-modules, Layer 1, Layer 2, and results after changing:

- Pando significance filtering;
- `top_k_neighbors`;
- `min_shared_tfs`;
- `min_tf_jaccard`;
- `max_targets_per_tf`;
- `expansion_mode`;
- subsystem or reaction cross-reference annotations;
- GEM GPR rules.

Stage 3 writes `merged_meta_modules.rds`. This object is a deduplicated biological reaction catalogue and contains no FASTCORE support.

### Rerun Stage 4 onward

Rerun Layer 1, Layer 2, and results after changing:

- `regulatory_alpha`;
- `tau`;
- gene half-saturation;
- metacell RNA or ATAC evidence while retaining the same Stage 3 catalogue.

### Rerun Stage 5 onward

Rerun only Layer 2 and results after changing:

- medium composition or exchange bounds;
- `target_direction`;
- LP solver or scoring `time_limit`;
- `omega`;
- global FASTCORE controls in `layer2_args$model_params`:
  - `completion_time_limit`;
  - `fastcore_epsilon`;
  - `max_support_reactions`;
  - `strict`.

## Rebuild only the medium-specific union GEMs

```r
step5_new <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = new_medium_scenarios,
  outdir = "RegCompass_restart/05_layer2",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 900,
    model_params = list(
      completion_time_limit = 900,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2500,
      strict = TRUE
    )
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

This reuses:

- the same biological meta-module definitions;
- the same merged complete-GPR target set;
- the same Layer 1 reaction-support matrix.

It reconstructs one union GEM for each new medium and performs one global FASTCORE completion per medium.

## Diagnose union-GEM completion

```r
summary <- step5_new$model_cache_summary
summary[, c(
  "medium_scenario",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "target_status",
  "file",
  "file_checksum"
)]

models <- lapply(summary$file, readRDS)

diagnostics <- do.call(
  rbind,
  Map(function(model, medium) {
    x <- model$closure_diagnostics
    x$medium_scenario <- medium
    x
  }, models, summary$medium_scenario)
)

table(diagnostics$medium_scenario, diagnostics$completion_status)
```

Interpretation:

- `already_feasible`: the merged biological catalogue already supported the target under that medium;
- `global_fastcore_completed`: global FASTCORE added supporting reactions;
- `parent_blocked`: the target direction is infeasible in the medium-constrained parent GEM;
- `unresolved`: completion failed under the requested constraints;
- `no_allowed_direction`: the original GEM bounds do not permit the requested direction.

## Compare support sets across media

```r
support_by_medium <- do.call(
  rbind,
  Map(function(model, medium) {
    ids <- model$reaction_meta$reaction_id[
      model$reaction_meta$global_fastcore_support %in% TRUE
    ]
    data.frame(
      medium_scenario = medium,
      reaction_id = ids,
      stringsAsFactors = FALSE
    )
  }, models, summary$medium_scenario)
)

table(support_by_medium$medium_scenario)
```

Differences in these support sets are structural consequences of different medium constraints. They should not be interpreted as condition-specific expression effects.
