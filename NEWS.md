# RegCompassR 2.4.12

- Removes absolute Pando `estimate` magnitude from the canonical Layer-1 regulatory penalty projection while retaining fitted estimates for inference diagnostics, significance testing, direction and audit.
- Uses `sign(estimate)` only after the existing fit-status, estimability and BH-adjusted-P gates in both condition-Pando and standard-Pando routes.
- Converts each paired-cell TF-by-ATAC predictor to a positive-scale-invariant bounded activity using a cell-type-wide robust scale before exact SuperCell membership aggregation, preventing arbitrary RNA/ATAC feature units from replacing the removed coefficient magnitude.
- Averages signed activities across each target/condition's active Pando edges so targets or conditions with more significant edges do not receive a larger regulatory score solely from network degree.
- Keeps the existing target-level robust calibration, bounded regulatory odds correction, RNA support, GPR aggregation, COMPASS reaction cost and Layer-2 LP unchanged.
- Adds regression tests for estimate-magnitude invariance, positive predictor-rescaling invariance, target active-edge-degree control, exact membership aggregation and sign-only routing, and updates the mathematical specification and tutorials.

# RegCompassR 2.4.7

- Removes finite structural time limits from the default CORDA2 route. CORDA2 now always runs with `completion_time_limit = Inf`, and supplying `model_params$completion_time_limit` with CORDA2 is rejected so a long Human-GEM reconstruction cannot be silently truncated by a persistent solver clock.
- Keeps `completion_time_limit` available for supplementary non-CORDA2 completion such as FASTCORE, and synchronizes the Layer 2 tutorials and function reference with this route-specific contract.
- Adds a regression test covering unlimited CORDA2 runtime and retained non-CORDA2 time-limit controls.
- Adds a post-analysis tutorial helper for same-cell-type Top-N reaction bar charts, defaulting to Top 20 by median support score, using formal reaction names or direction-specific formulas, and distinguishing RNA-only from active multiome-supported reactions with an any-metacell rule.
- Parallelizes independent FASTCC + FASTCORE structural reconstructions across cell-type-by-medium tasks through the active Layer 2 `BPPARAM` instead of running the cache builder as nested serial loops.
- Limits package-created structural worker pools to the number of independent FASTCORE tasks, enables dynamic task scheduling, and keeps solver libraries single-threaded inside each R worker to avoid nested oversubscription.
- Atomically checkpoints each completed FASTCORE model, returns only lightweight cache metadata to the controller, drops task-local model and table objects, runs full garbage collection in every worker, and releases package-managed worker processes after the structural batch.
- Reports requested/active structural workers, task count, dispatch policy, completion method, and worker-cleanup policy on the FASTCORE model cache.
- Reports Layer 2 as a 12-part end-to-end workflow instead of a single 0/1 progress item.
- Prints task-scoped progress for every cell-type-by-medium model, including COMPASS medium setup, FASTCC/FASTCORE phases or original CORDA2 Steps 1, 2.1, 2.2 and 3, target closure, directional vmax, primary penalty scoring and RNA-only control scoring.
- Persists merged `layer2_task_progress.tsv` and latest-state `layer2_task_status.tsv` diagnostics across controller and parallel worker processes, including the actual interrupted step on failure.
- Keeps numerical algorithms, model caches, parallel task decomposition and Layer 2 output schemas unchanged.
- Aligns medium application with the original COMPASS exchange-bound sequence: all exchange uptake directions receive the shared cap, omitted medium reactions remain capped rather than closed, and explicit medium rows override only the intended uptake direction for generated biological media.
- Detects uptake direction from exchange stoichiometry so both `metabolite -> boundary` and `boundary -> metabolite` conventions are handled without restricting the secretion direction.
- Keeps reaction-level custom `lb`/`ub` rows backward compatible, while `available = FALSE` closes uptake only.
- Expands the medium audit to allow COMPASS baseline changes on any annotated exchange reaction while continuing to reject bound changes on internal reactions.
- Versions the full-GEM medium fingerprint so existing caches created under the previous closed-unlisted semantics are not reused.

# RegCompassR 2.4.0

- Defines three mutually exclusive Layer 2 structural routes: add-only FASTCORE, original MATLAB CORDA2, and a COMPASS-style full GEM.
- Applies medium scenarios as exchange-reaction bound changes only; medium handling itself never removes reaction or metabolite columns in any of the three routes.
- Retains every reaction and requested target in full-GEM mode, computes directional maximum flux under the medium, and skips the penalty LP for medium-infeasible directions.
- Records `model_completion = "none"` for the full-GEM route and rejects FASTCORE/CORDA2 controls instead of silently accepting or mislabelling them.
- Fingerprints the exact medium bounds used by full-GEM caches so same-named media with different exchange limits cannot reuse stale models.
- Uses the latest default branches of Pando_regcompass and SuperCell_Seurat_V4 without fixed revisions.
- Routes each retained cell type independently to condition GRN or standard Pando according to its retained condition count.
- Uses the current Pando condition-fit API.
- Removes retired projection and penalty fields from Layer 1, Layer 2, documentation, and result schemas.
- Consolidates mathematical definitions in `docs/mathematical-model.md`.
- Reduces user tutorials and continuous integration to the current supported workflow.
