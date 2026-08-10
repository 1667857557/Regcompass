# CORDA2 performance notes

This page records runtime-facing rules only. The mathematical and structural invariants are maintained in [mathematical-model.md](mathematical-model.md).

- RegCompass uses the single top-level `workers` cap.
- Parallel work is restricted to independent tasks supported by the current CORDA2 stage; nested solver/thread pools are disabled.
- Each CORDA2 step is reduced before the next dependent step begins.
- CORDA2 structural reconstruction has no user-facing finite completion timeout.
- After the single CORDA2 reconstruction, the final GEM is scored directly; there is no second closure-LP phase to parallelize.

Use `step_progress.tsv`, `step_timing.tsv`, and the Layer 2 diagnostics in the returned object for runtime diagnosis.
