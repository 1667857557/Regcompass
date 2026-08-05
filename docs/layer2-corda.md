# Corrected Python CORDA2 reconstruction in Layer 2

`rc_regcompass_step_layer2()` keeps compact add-only FASTCORE as the default.
For `model_mode = "meta_module_gem"`, use:

```r
model_completion = "corda2"
```

The development aliases `"corda"` and `"corda_like"` are normalized to the
same CORDA2 implementation. No second CORDA reconstruction algorithm is loaded.

The implementation follows `resendislab/corda` at commit
`c02e06d50606bf93f23d8f2e6d6ade0e996ca70e`, including its five directional
confidence levels, redundant-path search and build stages. Three source-level
corrections are applied explicitly:

1. remaining medium-confidence directions are tested by **maximizing** their
   flux rather than minimizing a positive target coefficient;
2. the opposite copy of a reversible target is fixed to zero during target
   assessment, preventing a zero-net-flux forward/reverse self-cycle;
3. forward and reverse variables receive costs from their own current
   directional confidence rather than both inheriting the forward variable's
   current confidence after directional promotion.

The pinned Python package uses `np.in1d`, which is unavailable in NumPy 2.x.
Reference CI therefore runs that historical source with `numpy<2`; the R
implementation has no NumPy dependency.

## Recommended call

```r
step5_corda2 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_corda2",
  model_mode = "meta_module_gem",
  parallel = TRUE,
  BPPARAM = layer2_bp,
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      model_completion = "corda2",
      completion_time_limit = 3000,
      strict = TRUE,
      corda2_redundancies = 3L,
      corda2_penalty_factor = 100,
      corda2_support = 5L,
      corda2_cost_increase = 1.01,
      corda2_target_flux = 1,
      corda2_flux_tolerance = 1e-8,
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20,
      corda_include_evidence_outside_modules = TRUE,
      corda_max_medium_confidence_reactions = Inf
    )
  )
)
```

`fastcore_epsilon`, `max_support_reactions`, and `corda_seed` are not CORDA2
algorithm parameters. CORDA2 is deterministic and uses multiplicative cost
updates, not randomized cost perturbations.

## Parent-model contract

For each cell type and medium, CORDA2 receives the complete GEM after applying
the requested medium bounds. It does not run FASTCC first and does not block a
reaction solely from a generic demand, sink, exchange or artificial-support
role annotation.

As in the Python implementation, every direction that is open in the supplied
model is normalized to an internal upper bound of `1e6`. Closed directions stay
closed and positive directional lower bounds are preserved. The final retained
model restores the original reaction IDs, stoichiometry, annotations, GPR rules
and original bounds.

## RegCompass evidence mapping

Reaction evidence is calculated independently within each broad cell type:

\[
E_r=(1-w)\max\{Q_r^{RNA},Q_r^{multiome}\}+wA_r^{regulatory},
\]

where:

- `Q_RNA` is the within-cell-type percentile of median RNA-only reaction
  capacity;
- `Q_multiome` is the within-cell-type percentile of median multiome reaction
  capacity;
- `A_regulatory` is the median GPR regulatory-support fraction;
- the default regulatory weight is `w = 0.20`.

RegCompass maps reaction classes to Python CORDA2 confidence values:

| RegCompass class | Directional confidence | Meaning |
|---|---:|---|
| merged core reaction | `3` | high confidence |
| non-core module reaction | `2` | medium confidence |
| high-evidence outside-module reaction | `1` | low confidence |
| remaining low-evidence reaction | `-1` | absent / penalized |
| remaining reaction | `0` | unknown / free until final completion |

Forward and reverse copies inherit the same initial reaction confidence, then
are updated independently during reconstruction. A final reaction is retained
when at least one allowed directional copy reaches confidence `3`.

## Directional representation

Every allowed reaction direction is represented by a nonnegative variable:

\[
v_r=v_r^+-v_r^-.
\]

When direction `x` is assessed, its opposite directional copy is fixed to zero
and the target is forced to carry at least `corda2_target_flux`:

\[
v_x\ge t,\qquad v_{opposite}=0.
\]

The opposite-direction constraint is required because otherwise a reversible
reaction could satisfy the target through equal forward and reverse flux while
producing zero net stoichiometric flux.

## Associated-path LP

For one signed target, CORDA2 constructs costs from each directional variable's
current confidence:

\[
c_i=\begin{cases}
1,& q_i\in\{1,2\}\text{ and medium penalties are enabled},\\
P,& q_i=-1,\\
0,& q_i\in\{0,3\},
\end{cases}
\]

where the default absent-reaction penalty is `P = 100`.

The LP is:

\[
\begin{aligned}
\min_v\quad & \sum_i c_i v_i\\
\text{s.t.}\quad &S_{split}v=0,\\
&l\le v\le u,\\
&v_x\ge t,\\
&v_{opposite}=0.
\end{aligned}
\]

A low-, medium- or absent-confidence directional variable is associated with
the target when its solution flux exceeds `corda2_flux_tolerance`.

### Redundant pathways

For each target, the LP is solved repeatedly up to
`corda2_redundancies` times. After a solve, every newly observed penalized
variable has its coefficient multiplied by:

\[
CI=\texttt{corda2_cost_increase},
\]

with default `CI = 1.01`. The next solve is therefore encouraged to use an
alternative path. Repetition stops when no new associated variable is found or
the redundancy limit is reached.

Redundancy iterations for one target are necessarily serial because iteration
`k+1` depends on variables selected in iteration `k`. Different signed targets
remain independent and are parallelized.

## CORDA2 build stages

### Stage 1: high-confidence associations

All confidence-3 directional variables are assessed with low and medium
variables penalized by `1` and absent variables penalized by
`corda2_penalty_factor`. Every associated directional variable is promoted to
confidence `3`.

### Stage 2a: absent support for low/medium targets

Every remaining confidence-1 or confidence-2 direction is assessed with medium
penalties disabled. Associated absent directional variables are counted across
signed targets. An absent direction is promoted when its occurrence count is at
least `corda2_support`; the default is `5`, matching the Python package.

### Stage 2b: independent low/medium feasibility

All still-absent directions are blocked. Each remaining confidence-1 or
confidence-2 direction is tested independently by solving:

\[
\max v_x.
\]

The direction is promoted when the optimum is strictly greater than
`corda2_target_flux`.

The Python source sets a positive coefficient while retaining minimization in
this step, which drives the solution toward zero. RegCompass corrects the sign
and records this correction in the model and Layer-2 completion contracts.

### Stage 3: unknown/free reaction completion

All unpromoted low/medium directions are blocked. Unknown directions are changed
to absent confidence and become penalized candidates. Every confidence-3 target
is assessed once, with medium penalties disabled and redundancy search disabled.
Associated formerly unknown directions are promoted to confidence `3`.

## Native C++ solver acceleration

LP solving dominates runtime. HiGHS is already a native C++ solver, so the
implementation does not duplicate the optimizer in custom Rcpp code.

Each worker:

1. constructs one persistent HiGHS model for its target chunk;
2. updates the complete objective vector and target bounds;
3. repeatedly solves with the simplex method;
4. reuses the solver state and simplex basis;
5. falls back to the existing one-shot LP interface if the persistent API is
   unavailable or fails.

Custom C++ should be added only if profiling on Human-GEM shows that sparse
index preparation or result aggregation, rather than LP solving, is a material
runtime bottleneck.

## Parallel scheduling

Parallelism is adaptive:

- when `cell type × medium` model count is at least the worker count, models are
  constructed in parallel and each worker runs one model's CORDA2 stages;
- when model count is smaller, models run sequentially and independent signed
  targets use the complete worker pool;
- redundancy iterations for one target remain serial;
- confidence changes are applied only after every target in a build stage has
  completed;
- one `BPPARAM` pool is retained across the complete Layer-2 call;
- each worker restricts HiGHS/Gurobi and BLAS/OpenMP internals to one thread.

This avoids nested process pools and solver oversubscription.

## Pipeline contracts

- CORDA2 is available only with `model_mode = "meta_module_gem"`.
- One structural model is built per cell type and medium and shared across its
  matching conditions and metacells.
- The primary multiome calculation and RNA-only control reuse exactly the same
  model file and checksum.
- Only retained core-reaction directions that can reach
  `corda2_target_flux` enter downstream `vmax`, penalty and score matrices.
- `strict = TRUE` uses `corda2_target_flux`, not the much smaller numerical
  association tolerance, when deciding whether a parent-feasible core direction
  was lost.
- Low-, medium-, absent- and unknown-reaction reconstruction does not enlarge
  the score matrix.
- CORDA2 model files are isolated under:

```text
model_cache/meta_module_gem/corda2/
```

## Scope limitation: Python `met_prod`

The Python class can add temporary mock reactions through `met_prod`. The
current RegCompass Layer-2 API supplies reaction targets and medium constraints,
not metabolite-production tasks, so `met_prod` is not exposed in PR #254. The
implementation must not claim Python task parity until an explicit metabolite
or reaction-task input contract is added and tested.

## Audit output

Each CORDA2 model stores:

- the Python reference repository and commit;
- all initial and final directional confidence values;
- reaction-level inclusion stages;
- per-target redundancy counts;
- directional absent-support counts;
- impossible signed targets;
- solver runtime and persistent-solver fallback counts;
- all three intentional source corrections;
- compact task diagnostics;
- normalized association edges with both `associated_variable_id` and
  `associated_reaction_id`;
- original parent-model and downstream core-closure diagnostics.

Synthetic CI verifies the pinned Python behavior in a NumPy-compatible
environment and independently verifies R LP solutions, redundancy search,
support thresholds, independent medium promotion, final free-reaction
completion, reversible-target self-cycle prevention, serial/Multicore equality,
persistent/one-shot equality and installed-package SnowParam execution.
