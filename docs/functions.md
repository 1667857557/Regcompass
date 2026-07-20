# Public functions in RegCompassR 1.7.0

## `rc_prepare_gem()`

Prepare the supported human or mouse genome-scale metabolic model.

## `rc_prepare_human2_gem()` and `rc_prepare_mouse_gem()`

Prepare the pinned Human-GEM 2.0.0 or Mouse-GEM 1.8.0 model.

## `rc_make_medium_scenarios()`

Construct one shared extracellular medium. Condition-specific medium rows are
not accepted by shared-model scoring.

## `rc_run_regcompass()`

Run the canonical condition-pooled workflow.

Required design:

- `sample_col`: original biological-sample identifier; each sample maps to one condition;
- `condition_col`: condition used for pooling, ranking and optional comparison;
- `celltype_col`: cell type kept pure during pooling;
- `fragment_files = FALSE`;
- one or more samples per condition;
- one or more conditions.

A single condition produces `reaction_ranking` and no contrast. Multiple conditions
produce one ranking per condition and all pairwise descriptive contrasts within
each cell type. Pooled metacells are not treated as independent biological
replicates.

Main argument bundles:

- `metacell_args`: SuperCell2 parameters such as `gamma` and minimum stratum size;
- `pando_args`: `initiate_grn()`, motif and `infer_grn()` parameters;
- `layer1_args`: `regulatory_alpha`, `gene_half_saturation`, `tau`, and local FASTCORE options;
- `layer2_args`: solver, target direction, time limit and shared-model options.

Structural model selection:

- `model_mode = "meta_module_gem"`: score the shared union of locally completed condition-specific meta-modules;
- `model_mode = "full_gem"`: score the same targets and penalties in the shared full GEM.

Both modes use the same Layer 1 evidence model, medium and target-flux fraction.

Structural contract:

- core reactions require at least one complete GPR isozyme group;
- biological membership may expand through the core reaction's subsystem and
  shared KEGG, Reactome, or master-Rhea reaction identifiers;
- no reaction is added by metabolite sharing, stoichiometric adjacency, or a
  one-hop rule;
- there is no `include_one_hop` or metabolite-degree control in the API;
- local FASTCORE may add only reactions required for flux feasibility, reported
  separately from annotation-defined biological membership.

Evidence contract:

- Pando coefficients weight cell-type-referenced peak-accessibility deviations;
- the modifier updates bounded RNA support before GPR aggregation;
- protein-complex AND uses normalized Boltzmann soft-min;
- isozyme OR is additive;
- expression-linked reactions use `1 / (1 + log2(1 + E_multiome))`;
- only exchange, demand, sink and artificial-support reactions receive fixed
  structural costs.

## `rc_run_regcompass_one_shot()`

Prepare the species GEM and medium when omitted, then delegate to
`rc_run_regcompass()`. It uses the same condition-pooled evidence architecture,
selectable structural model modes and annotation-only meta-module expansion.
