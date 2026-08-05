# Exact Python CORDA2 semantics in Layer 2

`rc_regcompass_step_layer2()` keeps compact add-only FASTCORE as the default.
The optional CORDA2 route is selected only for
`model_mode = "meta_module_gem"`:

```r
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
    corda_medium_confidence_threshold = 0.75,
    corda_negative_confidence_threshold = 0.10,
    corda_regulatory_weight = 0.20,
    corda_include_evidence_outside_modules = TRUE,
    corda_max_medium_confidence_reactions = Inf
  )
)
```

The implementation follows `resendislab/corda` at commit
`c02e06d50606bf93f23d8f2e6d6ade0e996ca70e` **as written** for
`met_prod = NULL`. It does not substitute the 2016 paper algorithm and does not
apply the former PR-specific “corrections.” The compatibility names `"corda"`
and `"corda_like"` route to this same implementation.

## Exact-source scope

The Python constructor exposes three controls:

| Python argument | RegCompass argument | Default |
|---|---|---:|
| `n` | `corda2_redundancies` | `3` |
| `penalty_factor` | `corda2_penalty_factor` | `100` |
| `support` | `corda2_support` | `5` |

The following values are fixed in the source and are not configurable:

| Source symbol | Value | Role |
|---|---:|---|
| `UPPER` | `1e6` | normalized open directional bound |
| `CI` | `1.01` | multiplicative cost increase |
| `tflux` | `1` | target lower bound and iteration-2 comparison |

`corda2_cost_increase`, `corda2_target_flux`, `corda2_flux_tolerance`, and
`corda_seed` are rejected. CORDA2 has no randomized paper-noise process. The
association threshold is the active solver's feasibility tolerance, matching
`model.solver.configuration.tolerances.feasibility` in Python. RegCompass uses
`1e-7` for HiGHS/GLPK and `1e-6` for Gurobi.

The Python `met_prod` feature creates temporary mock metabolic reactions.
RegCompass Layer 2 currently supplies reaction targets and medium constraints,
not temporary metabolite-production reactions, so exactness is claimed only for
`met_prod = NULL`.

## Model initialization

For every reaction, both COBRA-style nonnegative variables are created:

\[
v_r = v_r^+ - v_r^- , \qquad v_r^+,v_r^-\ge 0.
\]

This includes closed directions. Reaction bounds are normalized in reaction
order exactly as in the Python constructor:

\[
l_r < -\tau \Rightarrow l_r=-10^6,
\qquad
u_r > \tau \Rightarrow u_r=10^6,
\]

where \(\tau\) is the active solver feasibility tolerance. Positive directional
lower bounds and directions whose bounds lie within \(\tau\) are preserved.
The final reconstructed GEM restores the original medium-constrained reaction
bounds.

The complete medium-constrained input GEM is copied before normalization. No
FASTCC pruning, role-based blocking, or global feasibility precheck is added.
`fastcore_epsilon` and `max_support_reactions` are irrelevant to CORDA2.

## Confidence representation

The Python algorithm expects reaction confidence in
`{-1, 0, 1, 2, 3}`. RegCompass derives those inputs within each broad cell type:

| RegCompass evidence class | CORDA2 confidence |
|---|---:|
| merged core reaction | `3` |
| non-core merged meta-module reaction | `2` |
| selected high-evidence outside-module reaction | `1` |
| finite low-evidence reaction | `-1` |
| remaining reaction | `0` |

Both direction variables initially receive the reaction confidence. Directional
confidence is then mutated by the Python build sequence. A final reaction is
included when `max(conf_forward, conf_reverse) == 3`.

The confidence mapping is RegCompass input preparation; the CORDA2 solver logic
begins after this mapping.

## `associated()` mathematical and code contract

Penalties are computed once before the target loop. The Python source looks up
the **forward-variable confidence** by reaction ID and assigns the same cost to
both direction variables:

\[
c_{r,+}=c_{r,-}=\begin{cases}
1, & q_{r,+}\in\{1,2\}\text{ and medium penalties are enabled},\\
P, & q_{r,+}=-1,\\
0, & \text{otherwise},
\end{cases}
\]

where \(P=\texttt{penalty\_factor}\). This remains true even when forward and
reverse confidence have diverged during the build.

Targets are processed serially in Python dictionary order. For target variable
\(x\):

1. If `ub(x) < tolerance`, append `x` to `impossible`, set its confidence to
   `-1`, and continue.
2. Otherwise save its bounds, set
   \(l_x=\max(1,l_x)\) and \(u_x=10^6\).
3. Do **not** block the opposite reversible variable.
4. Solve

\[
\begin{aligned}
\min_v\quad & \sum_i c_i v_i\\
\text{s.t.}\quad &S_{split}v=0,\\
&l\le v\le u.
\end{aligned}
\]

5. An active variable is associated when

\[
v_i>\tau,\qquad q_i\in\{-1,1,2\},\qquad i\ne x.
\]

6. For each newly observed associated variable that has an objective
   coefficient, multiply that coefficient by the fixed `CI = 1.01`.
7. Repeat while a new association is found and the iteration count is below
   `n`.
8. Restore only the target bounds.

The opposite reversible variable remains available. Consequently, the exact
source permits equal forward/reverse flux to satisfy a target with zero net
reaction flux. This behavior is preserved because the requirement is source
parity, not an algorithmic correction.

## Exact build sequence

### Iteration 1 — high-confidence dependencies

All directional variables with confidence `3` are passed to `associated()` with
medium penalties enabled. The union of returned directional variables is then
promoted to confidence `3`.

### Iteration 2a — absent support for low/medium targets

All remaining confidence-`1` or confidence-`2` variables are passed to
`associated()` with medium penalties disabled. Returned variables currently at
confidence `-1` are counted with Python `Counter` semantics. A directional
variable is promoted when its count is at least `support`.

Every remaining confidence-`-1` variable is then constrained by

\[
u_i=\max(0,l_i).
\]

### Iteration 2b — source positive-coefficient minimization

The Python source iterates serially over all solver variables. For each variable
whose current confidence is `1` or `2`, it sets a positive objective coefficient
and minimizes:

\[
\min_v v_i.
\]

The variable is promoted only when the optimum objective is strictly greater
than the fixed `tflux = 1`.

This is intentionally **not** replaced by flux maximization. In ordinary
zero-lower-bound models the optimum is commonly zero, so this source behavior
usually does not promote an otherwise optional medium-confidence variable.

### Iteration 3 — block remaining low/medium and test free variables

In confidence dictionary order:

- confidence `1` or `2`: set the directional upper bound directly to zero;
- confidence `0`: change confidence to `-1`.

Then all current confidence-`3` targets are passed to `associated()` once with
medium penalties disabled and redundancy detection disabled. Returned variables
are promoted to `3`.

Finally:

- `impossible` is sorted and deduplicated;
- redundancy records are retained only for final confidence-`3` variables not
  in `impossible`;
- a reaction is included when either directional confidence is `3`.

## Execution order and parallelism

Target-level parallelism is disabled because it would not reproduce the source's
single-model mutation and solver order. Each CORDA2 model uses one persistent
solver instance and processes targets and iteration-2 variables serially.

Parallelism is allowed only across independent `cell type × medium` model
instances. This does not change any CORDA2 state because each task creates its
own model, confidence dictionary, objective, bounds, and solver instance.

HiGHS is restricted to one thread within each worker. The persistent native
HiGHS interface updates objective coefficients and bounds on the same solver
model; the one-shot path is retained only as a runtime fallback.

## RegCompass post-reconstruction scoring

CORDA2 reconstruction and RegCompass scoring are separate contracts. After the
exact CORDA2 reaction subset is built, RegCompass evaluates core-reaction
forward/reverse capacities on the restored original bounds and includes only
core directions reaching `tflux = 1` in downstream `vmax`, penalty, and score
matrices.

The original Python build records impossible targets and completes rather than
throwing. Therefore `strict = TRUE` is retained for API compatibility but does
not convert a CORDA2 impossible/closure result into a reconstruction error.
Closure failures are stored in diagnostics and excluded from scoring.

## Cache and audit output

One structural model is built per cell type and medium and shared across matching
conditions and metacells. Primary multiome scoring and the RNA-only control reuse
the same model path and checksum. Files are stored under:

```text
model_cache/meta_module_gem/corda2/
```

Each model records:

- pinned Python repository and commit;
- constructor controls and fixed source constants;
- solver feasibility tolerance;
- initial, intermediate, and final directional confidence;
- impossible directional targets;
- redundancy counts;
- absent-support counts;
- source-order execution metadata;
- normalized association edges containing both directional variable ID and
  reaction ID;
- restored-bound core-closure diagnostics.

## Validation

The dedicated CI job installs the pinned Python revision under `numpy<2`, runs
an executable Python source oracle, and compares the R implementation on the
same nondegenerate networks. Checks cover:

- constants and solver tolerance;
- both direction variables and bound normalization;
- upstream redundant-path example;
- unblocked reversible self-cycle;
- forward-confidence penalty coupling;
- absent-support counting;
- positive-coefficient minimization with forced and optional medium variables;
- iteration-3 free-reaction completion;
- serial target order;
- persistent versus one-shot HiGHS equivalence on the test networks;
- installed-package/SnowParam namespace behavior.

Solver implementations may choose different members of a degenerate optimal
face. The algorithm, objective, constraints, mutation order, thresholds, and
confidence updates are source-identical; CI comparison uses networks with stable
optimal associations across the pinned Python GLPK and R HiGHS backends.
