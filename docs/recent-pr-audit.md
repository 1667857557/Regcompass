# Audit of recent RegCompassR pull requests

This audit uses the current `main` branch after PRs #166–#171 as the implementation baseline. Where an older PR description or review comment conflicts with a later merged change, the later merged source and public contract are authoritative.

## Audit method

The review compared:

1. every export in `NAMESPACE` with its current function formals and implementation;
2. internal argument forwarding, class checks, workflow parameters, GEM fingerprints, ordered units, saved artifacts, and restart behavior;
3. README, the public function index, Tutorials 1–5, workflow/stage documentation, and Rd usage blocks;
4. tests introduced by PRs #166–#171;
5. unresolved review threads on PRs #166, #167, and #169;
6. the successful static and installed-package checks on the final PR #171 head.

The audit distinguishes source-level contract verification from biological validation on a representative paired RNA+ATAC dataset. The latter remains an end-to-end analysis responsibility rather than a substitute for deterministic package tests.

## Exported function audit

| Public function | Current role | Audit result |
|---|---|---|
| `rc_prepare_gem()` | species-aware bundled/download GEM entry point | bundled and downloaded paths now both persist to `save_rds`; pinned human/mouse provenance retained |
| `rc_prepare_human2_gem()` | Human-GEM preparation | active species-specific convenience entry point |
| `rc_prepare_mouse_gem()` | Mouse-GEM preparation | active species-specific convenience entry point; mouse symbols retain source case |
| `rc_bundled_gem_manifest()` | installed GEM provenance | manifest fields and bundled checksums match the portable-execution contract |
| `rc_download_species_gem()` | explicit update/download path | retained for pinned release updates and asset rebuilding |
| `rc_make_medium_scenarios()` | shared medium construction | consistent across one-shot and stepwise tutorials |
| `rc_parallel_config()` | backend diagnostics | automatic SOCK/multicore resolution remains the documented policy |
| `rc_run_regcompass()` | canonical six-stage workflow | exposes two worker budgets; positional `species` compatibility is preserved ahead of the later `progress` argument |
| `rc_run_regcompass_one_shot()` | species-aware complete-run wrapper | prepares GEM/medium and forwards to the canonical runner |
| `rc_regcompass_step_grn()` | condition-by-cell-type Pando GRNs | classed checkpoint, GEM fingerprint, optional `BPPARAM`, committed-artifact timing |
| `rc_regcompass_step_metacells()` | condition-only label-guided SuperCell2 metacells | checkpoint reuse now requires an identical input/parameter cache contract |
| `rc_regcompass_step_meta_modules()` | complete-GPR cores, annotation expansion, local FASTCORE | matching GRN/metacell provenance and group coverage required |
| `rc_regcompass_step_layer1()` | integrated RNA+ATAC reaction expression | validates workflow/GEM/metacell provenance and ordered units |
| `rc_regcompass_step_layer2()` | shared-union/full-GEM directional LP scoring | validates Layer 1 provenance, solver availability, exact units, and persistent model files |
| `rc_regcompass_step_results()` | annotated final result assembly | rejects legacy metacell artifacts before assigning current construction provenance |
| `rc_regcompass_step_target_union()` | direct database-linked non-core second-pass scoring | validates each selector independently and checks candidates against actual cached union models |
| `rc_build_reaction_annotations()` | reaction annotation and evidence construction | uses normalized GEM bounds, inferred roles, source-case genes, and conservative reaction-level evidence classes |
| `rc_attach_reaction_annotations()` | attach annotations to an existing result | older results cannot claim `RNA+ATAC` when reaction-capacity reconstruction is unavailable |
| `rc_test_condition_reactions()` | same-target condition comparisons | retains annotation/evidence when called from plotting helpers |
| `rc_select_gene_reactions()` | case-insensitive GPR-aware selection | returns the source symbol case rather than forcing mouse genes to uppercase |
| `rc_plot_condition_reaction()` | one-reaction multi-condition plot | passes the full annotated result into the statistics layer |
| `rc_plot_condition_gene_reactions()` | gene-associated plot collection | applies requested conditions to evidence selection and exposes/forwards `min_units` |

## PR claim verification

### PR #166: reaction annotation, evidence, statistics, and plots

The broad functionality was present, but nine unresolved review findings remained valid in current `main`:

- single-reaction plots discarded annotation context by passing only `microcompass`;
- gene-plot evidence filtering ignored the requested conditions;
- `rc_build_reaction_annotations()` lacked its own roxygen export declaration;
- unavailable reaction-capacity comparisons could fall back to gene changes and incorrectly claim `RNA+ATAC`;
- normalized bounds returned by `rc_validate_gem()` were discarded;
- mouse symbols were forced to uppercase;
- gene plots silently hard-coded `min_units = 5`;
- omnibus rows with no matched evidence defaulted to `structural/no-GPR`;
- reactions without an explicit role defaulted to `internal` instead of using the existing role annotator.

All nine are corrected in this follow-up.

### PR #167: target-union restart semantics

Two review findings had already been fixed by later strict stage contracts: cross-run stage mixing and loss of blocked-direction diagnostics. Two remained valid:

- target availability was inferred from a legacy global membership table rather
  than the exact cached union GEM files;
- an invalid `core_genes` selector could be ignored when a valid explicit reaction ID was also supplied.

The second-pass definition now reads each cached medium-specific GEM, uses their reaction intersection as the reusable target universe, and validates gene and reaction selectors independently.

### PR #168: strict stage contracts and direct non-core scoring

Implemented in current source:

- formal Layer 1 and Layer 2 classes;
- workflow parameter, GEM fingerprint, core-set, and ordered-unit validation;
- direct KEGG/Reactome/master-Rhea mapping;
- exclusion of same-subsystem-only, recursive/transitive, generic-neighbour, and previously scored global-core targets;
- exact reuse of medium-specific global union-GEM files.

This follow-up strengthens the exact-reuse claim by making the cached files, rather than the pre-completion membership table, authoritative for target availability.

### PR #169: portable execution, bundled GEMs, progress, and timing

Two review findings were already resolved by later manifest/schema changes. Three remained valid:

- bundled models were returned from `inst/extdata` without being copied to a requested `save_rds` path;
- adding `progress` before `species` changed positional calls to `rc_run_regcompass()`;
- a stage could write `status = success` before its final RDS/export completed.

Bundled and downloaded models now share the same persist-and-revalidate path; `species` again precedes `progress`; known public stages finalize success only after the expected final RDS is newly committed.

### PR #170: layered parallelism and worker lifecycle

Implemented in current source:

- defaults of six upstream workers and thirty Layer 2 workers;
- automatic SOCK workers on Windows and multicore workers on Linux/macOS;
- serial behavior for one worker;
- stage-scoped create/start/stop/release lifecycle;
- one numerical thread per outer task;
- Pando inner parallelism forced off in the canonical runner.

No current unresolved review finding required an algorithmic change.

### PR #171: tutorials and `gamma = 30`

Implemented in current source and documentation:

- `gamma = 30` in canonical and low-level defaults;
- a one-shot tutorial distinct from a true six-stage stepwise tutorial;
- a targeted gene/reaction remapping tutorial;
- a condition differential-analysis tutorial.

The stale function-index label was corrected, and Tutorials 3–5 were synchronized with cache contracts, actual cached-model target availability, conservative evidence resolution, and explicit `min_units` forwarding.

## Additional problems found by the function-by-function audit

### Stage 2 checkpoints were reused by file existence only

The low-level builder returned completed stratum directories whenever expected files existed and `overwrite = FALSE`. A changed `gamma`, label column, cell set, condition assignment, assay content, or dimensional setting could silently reuse stale memberships.

The canonical stage now writes `condition_metacell_cache_contract.rds`. It records ordered cells and metadata, complete RNA/ATAC assay-content projections computed in O(nnz) without dense materialization, the exact PCA/LSI embedding values used by SuperCell2, the construction label, reductions/dimensions, `gamma`, seed, and metacell thresholds. Existing checkpoints are reused only when the contract is identical. Legacy or changed caches require `metacell_args = list(overwrite = TRUE)`.

### Legacy metacell objects could receive current provenance

A manually loaded older classed object could pass a class-only check and later be described as condition-only label-guided output. Downstream stages now require the current input-design record and cache-contract schema before current provenance is assigned.

### Dependency declarations and the exact runtime gate were conflated

`DESCRIPTION` uses the legal minimum-version syntax supported by R: SeuratObject `>= 4.1.4`, Seurat `>= 4.4.0`, and Signac `>= 1.11.0`. R package metadata does not support an equality operator in dependency fields. The exact validated stack remains recorded in the three `Config/RegCompassR/*Version` fields and enforced by `.onLoad`, which rejects versions other than 4.1.4, 4.4.0, and 1.11.0. Tests now verify both layers instead of writing an invalid `DESCRIPTION`.

### Public documentation lagged the implementation

The function index, generated Rd files, targeted-remapping tutorial, and condition-analysis tutorial were synchronized with the actual signatures and failure policies.

## Deliberately unchanged

This follow-up does not change:

- Pando fitting equations or thresholds;
- SuperCell2 clustering mathematics or the `gamma = 30` default;
- complete-GPR aggregation;
- biological meta-module expansion;
- FASTCORE feasibility completion;
- reaction-expression or penalty equations;
- directional LP objectives or solver settings;
- direct database-link definition for target union;
- Wilcoxon/Kruskal-Wallis statistics;
- worker counts or backend policy.

## Regression coverage

The added or expanded tests verify that:

- missing or changed metacell cache contracts fail closed;
- legacy metacell stages are rejected downstream;
- legal dependency minimums and exact Config pins match the runtime gate;
- bundled GEMs are persisted to requested paths;
- `species` remains positionally ahead of `progress`;
- known stages report success only after their final RDS exists;
- unavailable reaction-capacity evidence cannot claim `RNA+ATAC`;
- missing omnibus evidence is `unknown/unavailable`;
- normalized bounds, inferred roles, and mouse symbol case are preserved;
- condition plots retain reaction annotations;
- gene plots apply requested conditions and `min_units`;
- each target-union selector is validated independently;
- reactions present in cached union models but absent from the membership table remain eligible when directly database-linked.
