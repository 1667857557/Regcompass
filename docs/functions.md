# Public functions in RegCompassR 1.7.0

## `rc_prepare_gem()`

Prepare the supported human or mouse genome-scale metabolic model.

## `rc_prepare_human2_gem()` and `rc_prepare_mouse_gem()`

Prepare the pinned Human-GEM 2.0.0 or Mouse-GEM 1.8.0 model.

## `rc_make_medium_scenarios()`

Construct one shared extracellular medium. Condition-specific medium rows are
not accepted by shared-GEM scoring.

## `rc_run_regcompass()`

Run the canonical condition-pooled workflow.

Required design:

- `sample_col`: original biological-sample identifier;
- `condition_col`: condition used for pooling and comparison;
- `celltype_col`: cell type kept pure during pooling;
- `fragment_files = FALSE`;
- `inference_unit = "metacell"`.

Main argument bundles:

- `metacell_args`: SuperCell2 parameters such as `gamma` and minimum stratum size;
- `pando_args`: `initiate_grn()`, motif and `infer_grn()` parameters;
- `layer1_args`: `regulatory_alpha`, `gene_half_saturation`, `tau`, and local FASTCORE options;
- `layer2_args`: solver, target direction, time limit and shared-model options.

Structural contract:

- core reactions require at least one complete GPR isozyme group;
- biological membership may expand through the core reaction's subsystem and
  shared KEGG, Reactome, or master-Rhea reaction identifiers;
- no reaction is added by metabolite sharing, stoichiometric adjacency, or a
  one-hop rule;
- there is no `include_one_hop` or metabolite-degree control in the API;
- local FASTCORE may add only the reactions required for flux feasibility, and
  these are reported separately from annotation-defined biological membership.

## `rc_run_regcompass_one_shot()`

Prepare the species GEM and medium when omitted, then delegate to
`rc_run_regcompass()`. It uses the same annotation-only meta-module expansion and
has no metabolite-neighbour expansion interface.
