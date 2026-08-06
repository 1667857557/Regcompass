# RegCompassR 2.4.7

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
