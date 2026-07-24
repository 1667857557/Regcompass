# Audit of recent RegCompassR pull requests

This audit uses the current `main` branch after PRs #168–#171 as the implementation baseline. Where older PR descriptions conflict with later changes, the later merged behavior is treated as authoritative.

## Audit method

The review compared:

1. the exported API in `NAMESPACE`;
2. public function formals and internal forwarding paths;
3. stage classes, workflow parameters, GEM fingerprints, ordered scoring units, and saved artifacts;
4. README, the public function index, all five tutorials, the workflow vignette, and Rd contracts;
5. tests added for PRs #168–#171 and the successful checks on the final PR #171 head.

The audit distinguishes source-level contract verification from full biological runtime validation. The latter remains the responsibility of the installed-package GitHub Actions workflow and representative end-to-end datasets.

## Exported function audit

| Public function | Current role | Audit result |
|---|---|---|
| `rc_prepare_gem()` | species-aware bundled/download GEM entry point | consistent with bundled Human-GEM 2.0.0 and Mouse-GEM 1.8.0 documentation |
| `rc_prepare_human2_gem()` | Human-GEM preparation | retained as an active species-specific preparation API |
| `rc_prepare_mouse_gem()` | Mouse-GEM preparation | retained as an active species-specific preparation API |
| `rc_bundled_gem_manifest()` | installed GEM provenance | consistent with the portable-execution documentation |
| `rc_download_species_gem()` | explicit update/download path | retained as the documented low-level update path |
| `rc_make_medium_scenarios()` | shared medium construction | consistent with one-shot and stepwise tutorials |
| `rc_parallel_config()` | backend diagnostics | consistent with automatic SOCK/multicore resolution |
| `rc_run_regcompass()` | canonical six-stage workflow | exposes only `upstream_workers` and `layer2_workers`; stage-scoped pools and single-thread child policy are implemented |
| `rc_run_regcompass_one_shot()` | species-aware complete-run wrapper | prepares GEM/medium and forwards to the canonical runner |
| `rc_regcompass_step_grn()` | condition-by-cell-type Pando GRNs | public stage, classed checkpoint, GEM fingerprint, optional `BPPARAM` |
| `rc_regcompass_step_metacells()` | condition-only label-guided SuperCell2 metacells | public stage; this follow-up adds checkpoint-content and construction-provenance validation |
| `rc_regcompass_step_meta_modules()` | complete-GPR cores, annotation expansion, local FASTCORE | requires matching GRN/metacell stages and validates bidirectional group coverage |
| `rc_regcompass_step_layer1()` | integrated RNA+ATAC reaction expression | validates stage settings, GEM identity, and ordered metacell units |
| `rc_regcompass_step_layer2()` | shared-union/full-GEM directional LP scoring | validates Layer 1 provenance, solver availability, and exact unit order |
| `rc_regcompass_step_results()` | annotated final result assembly | validates all upstream stages; legacy metacell objects are now rejected before current provenance is assigned |
| `rc_regcompass_step_target_union()` | direct database-linked non-core second-pass scoring | matches the latest direct KEGG/Reactome/master-Rhea, non-core-only contract and reuses exact cached union GEMs |
| `rc_build_reaction_annotations()` | reaction annotation construction | consistent with Stage 6 annotation outputs |
| `rc_attach_reaction_annotations()` | attach annotations to existing results | consistent with the public interpretation API |
| `rc_test_condition_reactions()` | same-target condition comparisons | documented as metacell-level within-dataset inference |
| `rc_select_gene_reactions()` | GPR-aware scored-reaction selection | consistent with the targeted interpretation tutorial |
| `rc_plot_condition_reaction()` | one-reaction multi-condition plot | consistent with the differential-analysis tutorial |
| `rc_plot_condition_gene_reactions()` | gene-associated plot collection | consistent with the reaction annotation/statistics contract |

## Recent PR claim verification

### PR #168: strict stage contracts and direct non-core scoring

Implemented in current source:

- formal Layer 1 and Layer 2 classes;
- workflow parameter, GEM fingerprint, and ordered-unit validation;
- direct KEGG/Reactome/master-Rhea mapping;
- exclusion of same-subsystem-only, recursive, FASTCORE-only, generic-union, and previously scored global-core targets;
- exact reuse of medium-specific global union-GEM files.

### PR #169: portable execution, bundled GEMs, progress, and timing

Implemented in current source and documentation:

- bundled pinned human and mouse GEM assets;
- species-aware preparation paths;
- public-stage progress and `step_timing.tsv`;
- complete-run `00_execution_timing.tsv` and result timing metadata.

### PR #170: layered parallelism and worker lifecycle

Implemented in current source:

- defaults of six upstream workers and thirty Layer 2 workers;
- automatic SOCK workers on Windows and multicore workers on Linux/macOS;
- serial behavior for one worker;
- stage-scoped create/start/stop/release lifecycle;
- one internal numerical thread per outer task;
- Pando inner parallelism forced off in the canonical runner.

### PR #171: tutorials and `gamma = 30`

Implemented in current source and main documentation:

- `gamma = 30` in canonical and low-level defaults;
- a one-shot tutorial distinct from a true six-stage stepwise tutorial;
- a targeted gene/reaction remapping tutorial;
- a condition differential-analysis tutorial.

One stale public-index label still called Level 2 a saved-stage audit; this follow-up changes it to the true stepwise workflow and lists Levels 4 and 5 explicitly.

## Problems found and fixed in this follow-up

### 1. Existing metacell checkpoints were reused by file existence only

The low-level builder returned completed stratum directories whenever expected files existed and `overwrite = FALSE`. It wrote `run_params.yaml` after construction but did not compare it before reuse. A changed `gamma`, cell-type label column, cell set, condition assignment, assay content, or dimensional setting could therefore silently reuse stale memberships.

The canonical condition-metacell stage now writes `condition_metacell_cache_contract.rds` before construction. The contract covers ordered cells and metadata, RNA/ATAC assay fingerprints, the construction label, reductions/dimensions, `gamma`, seed, and metacell size thresholds. Existing checkpoints are reused only when the contract is identical. Legacy or changed caches require `metacell_args = list(overwrite = TRUE)`.

### 2. Legacy metacell artifacts could receive current provenance

Final-result assembly described the metacell policy as label-guided even when a manually loaded older classed object lacked evidence that it was constructed that way. The stage-class check previously verified only the class name.

Current downstream stages now require the condition-only, label-guided input-design record and the new cache-contract schema. Legacy or incompatible artifacts stop with an explicit rebuild instruction before Layer 1, meta-module construction, or Stage 6 can assign current provenance.

### 3. Declared dependency bounds contradicted the runtime gate

`DESCRIPTION` allowed newer SeuratObject, Seurat, and Signac versions with `>=`, while `.onLoad` rejected every version except 4.1.4, 4.4.0, and 1.11.0. The declared Imports now use exact equality constraints, matching installation instructions, package configuration, and runtime enforcement.

### 4. Public tutorial index wording was stale

`docs/functions.md` still labelled Level 2 as a saved-stage audit after PR #171 restored direct execution of all six stages. The index now calls it the true stepwise workflow and explicitly links the Level 4 and Level 5 tutorials.

## Deliberately unchanged

This follow-up does not change:

- GRN fitting equations or Pando thresholds;
- SuperCell2 clustering mathematics or the `gamma = 30` default;
- complete-GPR aggregation;
- meta-module biological expansion;
- FASTCORE feasibility completion;
- reaction-expression or penalty equations;
- directional LP objectives;
- target-union direct-link semantics;
- condition statistical tests;
- worker counts or backend policy.

## Added regression coverage

The follow-up tests verify that:

- checkpoints without a cache contract cannot be silently reused;
- changed contracts are rejected unless overwrite is explicit;
- legacy metacell stage objects are rejected by downstream contracts;
- exact dependency declarations match the runtime version gate;
- the public function index retains the true stepwise tutorial wording.
