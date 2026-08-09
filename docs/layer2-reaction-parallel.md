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

Condition-GRN and standard-Pando candidate screening are kept separate from the
post-fit penalty-entry decision.

For a common-dictionary condition GRN, `tf_cor` and `peak_cor` act during
candidate discovery. After the same frozen dictionary has been fitted in every
condition, a TF-peak-gene coefficient enters regulatory penalty projection only
when:

- the coefficient is estimable;
- BH-adjusted `padj < 0.05`;
- the fitted estimate is finite.

There is no second post-fit `abs(corr)` or raw `abs(estimate)` cutoff. The stored
`penalty_effect` is the fitted condition coefficient for eligible edges and zero
otherwise. The retained correlation/estimate threshold provenance constants are
zero and indicate that no additional post-fit effect-size gate is applied.

For standard Pando, the requested `tf_cor` is a candidate-screening floor and may
be raised by the sample-size-aware Pearson critical-correlation floor. Its
post-fit active-edge filter likewise uses adjusted `P < 0.05` and estimability
when available, rather than a second fixed raw effect-size threshold.

Pando compatibility is checked from required exported APIs and data contracts,
not from a package-version floor.
