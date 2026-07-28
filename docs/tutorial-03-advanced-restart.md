# Tutorial Level 3: restart, sensitivity, and diagnostics

RegCompass saves classed stage objects so downstream stages can be rerun without repeating unchanged work. Objects from incompatible runs are rejected through stored workflow and GEM provenance.

This tutorial uses the current six-stage API documented in
[functions.md](functions.md). Restart from the earliest stage whose persisted
contract changes.

## Load a completed stepwise run

```r
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_steps/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

## Restart boundaries

### Rerun Stage 1 onward

Rerun Pando and all downstream stages after changing:

- `species`;
- `min_abs_estimate` or `min_model_rsq`;
- `reference_condition`, `tf_cor`, `peak_cor`, alpha, lambda-path/CV settings,
  or other supported `pando_infer_args`;
- motif matrices or genome;
- `pando_initiate_args$regions`;
- RNA/ATAC matrices, condition metadata, or cell-type metadata.

When `pfm` is omitted, the canonical motif collection is Pando's `motifs` data object. Supplying a different `pfm` changes the fitted regulatory evidence and therefore requires Stage 1 onward to be rerun.

Changing the edge dictionary, eligibility mask, pooled TF×ATAC transform,
target scale, or selected lambda invalidates every downstream condition effect
and requires Stage 1 onward to be rerun.

The species-specific default region policies are:

```text
human = phastConsElements20Mammals.UCSC.hg38 ∪ SCREEN.ccRE.UCSC.hg38
mouse = phastConsElements20Mammals.UCSC.hg38 only
```

Changing `species`, overriding the region object, or changing either default region source requires Stage 1 onward to be rerun.

### Rerun Stage 2 onward

Rerun metacells and all downstream stages after changing cell membership, RNA/ATAC matrices, condition or cell-type metadata, `gamma`, or the RNA/ATAC reductions used for metacell construction.

### Rerun Stage 3 onward

Rerun reaction meta-modules and downstream stages after changing:

- a custom `meta_module_args$subsystem_table`;
- subsystem, KEGG, Reactome, or master-Rhea annotations;
- GEM GPR rules.

Pando filtering is configured in Stage 1. Stage 3 uses the resulting supported target-gene set and accepts an optional custom subsystem table.

Stage 3 always performs one ordered pass:

```text
core subsystem
→ KEGG/Reactome equivalence
→ master-Rhea equivalence
```

### Rerun Stage 4 onward

Rerun Layer 1 and downstream stages after changing:

- `regulatory_alpha`;
- `gpr_and_method` among `"min"`, `"median"`, and `"mean"`;
- gene half-saturation;
- metacell RNA or ATAC evidence.

The default GPR-AND method is `"min"`.

### Rerun Stage 5 onward

Rerun Layer 2 and results after changing:

- medium composition or exchange bounds;
- `target_direction`, `omega`, solver, or `flux_threshold`;
- `completion_time_limit`, `fastcore_epsilon`, `max_support_reactions`, or `strict` in `layer2_args$model_params`.

`completion_time_limit` controls only the FASTCORE/FASTCC work that constructs the medium-specific union GEM. Scoring LPs have no `time_limit` API and run after the union GEM has been completed and cached.

The complete preset list and custom-medium format are documented in [medium presets](medium-presets.md).

## Rerun Stage 4 with another COMPASS GPR-AND rule

```r
step4_mean <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_restart/04_layer1_mean",
  regulatory_alpha = 1,
  gpr_and_method = "mean",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

Use `"median"` or `"mean"` for sensitivity analysis. The canonical default remains `"min"`.

## Rebuild Stage 5

```r
new_medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "low_glucose",
  species = "human"
)

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

This reuses the Stage 3 reaction targets and Stage 4 support matrix. The new `completion_time_limit` affects union-GEM reconstruction only; the subsequent scoring phase is not time limited.

## Diagnose model completion

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

model <- readRDS(summary$file[[1]])
diagnostics <- model$closure_diagnostics
table(diagnostics$completion_status)
```

Common statuses are:

- `already_feasible`: the initial reaction set supports the target;
- `global_fastcore_completed`: FASTCORE added supporting reactions;
- `parent_blocked`: the target direction is infeasible in the medium-constrained parent GEM;
- `unresolved`: completion did not succeed under the requested construction limit;
- `no_allowed_direction`: the original GEM bounds block the requested direction.
