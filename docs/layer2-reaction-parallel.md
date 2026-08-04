# Stage 5 reaction-level parallel execution

Stage 5 uses task-level process parallelism. Every worker process keeps the
numerical libraries and LP solver at one internal thread; Gurobi is explicitly
called with `Threads = 1`.

## FASTCORE construction

Union-GEM construction uses a two-level, non-nested plan:

- when several cell-type/medium models are requested, independent
  `cell type × medium` FASTCORE builds are dispatched across workers;
- when only one model is requested, the parent, initial and final directional
  feasibility checks are dispatched as independent directional-reaction tasks;
- an outer FASTCORE worker forces its inner feasibility checks to serial mode,
  preventing nested worker pools.

Each completed union GEM is written atomically to the model cache. The worker
then drops the model object and runs full garbage collection before accepting
another task.

## Directional vmax

The shared structural `vmax` cache uses one task for each exact
`cell type × medium × reaction × direction` row. A task loads its union GEM,
solves one directional LP, atomically writes a reaction checkpoint, releases the
model and invokes full garbage collection. Completed checkpoints are reusable
within the same model cache.

## Step-2 penalty LP

The Step-2 task unit is one exact directional reaction. The worker loads the
matching union GEM once and evaluates that reaction across all matching
metacells. This avoids reloading the model for every metacell while allowing
multiple reactions to run concurrently.

After one reaction is complete, the worker atomically writes a small checkpoint,
removes the union GEM and temporary solver objects, and runs full garbage
collection. The controller reads each checkpoint into the final matrices and
removes the temporary file. Worker processes are reused rather than killed and
recreated; only reaction-specific memory is discarded between tasks.

## Stage 1 penalty-entry contract

A TF-peak-gene edge enters regulatory penalty projection only when all conditions
below hold:

- the coefficient is estimable;
- BH-adjusted `padj < 0.05`;
- `abs(corr) >= 0.1`;
- `abs(estimate) >= 0.01`.

For standard Pando, `corr` is the coefficient-table TF-target correlation. For a
common-dictionary condition fit that does not expose a coefficient-level `corr`,
RegCompass uses the dictionary's audited `max_abs_tf_target_cor` and records the
source in `corr_source`.

Pando compatibility is checked from required exported APIs and data contracts,
not from a package-version floor.
