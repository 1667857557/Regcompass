# Original CORDA reconstruction in Layer 2

`rc_regcompass_step_layer2()` keeps compact add-only FASTCORE as the default.
For `model_mode = "meta_module_gem"`, `model_completion = "corda"` runs the
reaction-reconstruction algorithm described by Schultz and Qutub (2016). The
former draft value `"corda_like"` is accepted only as a compatibility alias and
is normalized to `"corda"`.

The implementation contract is based on the paper's Materials and Methods
pseudocode and the author-distributed MATLAB function referenced as S2 File. If
a RegCompass engineering choice is not part of the original algorithm, it is
identified explicitly below.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_corda",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      model_completion = "corda",
      completion_time_limit = 3000,
      strict = TRUE,
      corda_gamma = 1e5,
      corda_kappa = 1e-2,
      corda_epsilon = 1,
      corda_n = 5L,
      corda_p = 2L,
      corda_seed = 1L,
      corda_flux_tolerance = 1e-8,
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20,
      corda_include_evidence_outside_modules = TRUE,
      corda_max_medium_confidence_reactions = Inf
    )
  )
)
```

`fastcore_epsilon` and `max_support_reactions` are FASTCORE controls. They do not
change CORDA reconstruction and should normally be omitted from a CORDA call.

## Input model contract

For each cell type and medium, CORDA receives the complete GEM after application
of the requested medium bounds. No FASTCC reaction deletion is performed and no
reaction is blocked solely because it is annotated as demand, sink, exchange or
artificial support. This follows the original function's contract: CORDA acts on
the model supplied to it, and every supplied reaction is assigned to HC, MC, NC
or OT.

The medium constraint is a RegCompass input transformation performed before
CORDA. It is not presented as part of the original CORDA algorithm.

## Confidence mapping

CORDA accepts reaction confidence classes and leaves their definition to the
caller. RegCompass maps evidence within each broad cell type.

RNA-only and multiome reaction capacities are summarized across matching
metacells and converted to within-cell-type percentile ranks. The RegCompass
evidence score is

\[
E_r=(1-w)\max\{Q_r^{RNA},Q_r^{multiome}\}+wA_r^{regulatory},
\]

where `A_regulatory` is the median GPR regulatory-support fraction and the
default `w` is `0.20`.

The initial classes are:

- **HC:** merged core reactions;
- **MC_module:** non-core reactions in the merged cell-type meta-module;
- **MC_evidence:** optional reactions outside the module with evidence at or
  above `corda_medium_confidence_threshold`;
- **NC:** remaining reactions with finite evidence at or below
  `corda_negative_confidence_threshold`;
- **OT:** every remaining input-model reaction, including reactions without
  available omics evidence.

MC reactions are flexible and are not automatically retained.

## Exact directional transformation

The original dependency code does not split the reaction currently being
tested. It constrains that original variable to `+epsilon` or `-epsilon` and
splits the other reversible reactions.

RegCompass pre-splits the model once for solver reuse. To remain mathematically
equivalent, when a forward target is tested its reverse copy is fixed to zero;
when a reverse target is tested its forward copy is fixed to zero. Without this
constraint, a reversible reaction could satisfy the target by sending equal
forward and reverse flux through itself, producing zero net stoichiometric flux.

Original reaction bounds determine which directional copies exist. The final
model is converted back to original reaction IDs, stoichiometry, bounds,
annotations and GPR rules.

## Dependency-assessment LP

Let `Y` be the undesirable reaction group for the current stage. After splitting
non-target reversible reactions, every directional variable is nonnegative. The
original code adds one pseudo-metabolite whose production coefficient is the
reaction cost, adds a cost-consuming reaction, and minimizes that sink flux.
Eliminating the pseudo-metabolite balance gives the equivalent LP used here:

\[
\begin{aligned}
\min_v\quad &\sum_i c_i v_i\\
\text{s.t.}\quad&S_{split}v=0,\\
&l\le v\le u,\\
&v_x\ge\epsilon,
\end{aligned}
\]

where the opposite target direction is fixed to zero. For every directional
reaction,

\[
c_i=\begin{cases}
\gamma+U_i(0,\kappa), & i\in Y,\\
U_i(0,\kappa), & i\notin Y.
\end{cases}
\]

Thus CORDA minimizes the **combined flux through costly reactions**, not the
number of costly reactions. A reaction in `Y` is associated with the target when
its selected directional flux is nonzero. `corda_flux_tolerance` supplies the
numerical nonzero threshold; `corda_epsilon` is kept separate and defines the
forced target magnitude.

The original algorithm draws fresh uniform noise for every dependency
assessment. RegCompass preserves that distribution but derives each task's
random stream deterministically from `corda_seed`, stage, target direction and
repeat. This is an engineering adaptation required so serial, Multicore and SOCK
execution return the same task-level random draws regardless of scheduling. It
is not claimed to reproduce a particular unseeded MATLAB random-number stream
bit for bit.

The paper defaults are:

- `corda_gamma = 1e5`;
- `corda_kappa = 1e-2`;
- `corda_epsilon = 1`;
- `corda_n = 5` dependency assessments per allowed target direction;
- `corda_p = 2` distinct MC reactions required to promote a shared NC reaction.

## Three reconstruction stages

### Stage 1: HC associations

All HC reactions begin in RE. Every allowed HC direction is assessed `n` times.
MC directions have base cost `sqrt(gamma)`, NC directions have base cost
`gamma`, and HC/OT directions have base cost zero. Any MC or NC reaction
associated with an HC target in any repeat is moved to RE.

### Stage 2: flexible MC/NC core

Each remaining MC direction is assessed `n` times with NC as the costly group.
An NC reaction is moved to RE when it is associated with at least `p` **distinct
MC reactions**. All remaining NC reactions are then blocked. Every remaining MC
reaction is tested in every allowed direction and is moved to RE if at least one
direction can carry `corda_epsilon`.

The algorithm records excluded MC reactions and their Stage-2 NC associations
for subsequent curation.

### Stage 3: OT associations

All remaining MC and NC reactions are blocked. Every allowed direction of every
RE reaction is assessed `n` times with OT as the costly group. Every associated
OT reaction is moved to RE. The resulting RE set defines the final reconstruction.

Confidence updates occur only after all tasks in a stage finish. This stage
barrier preserves the original algorithm's reaction-order independence.

## Optional metabolic tasks

The original MATLAB workflow can temporarily add metabolite sinks as HC tests,
removing each test before another reaction or task is assessed. The current
RegCompass Layer-2 input contract supplies reaction targets rather than temporary
metabolite-task definitions, so this optional original-code feature is not yet
exposed. The implementation must not claim metabolic-task preservation unless a
future explicit task input and isolated temporary-reaction contract are added.

## Native C++ acceleration and parallelism

Repeated LP solution dominates runtime. HiGHS already provides the native C++
solver, so RegCompass does not duplicate FBA in custom Rcpp code. Each worker:

1. constructs one persistent HiGHS model for its task block;
2. updates only objective coefficients and target bounds;
3. re-solves with simplex-basis reuse;
4. falls back to the existing one-shot LP path if the persistent API is absent
   or fails.

When the number of cell-type-by-medium models is at least the worker count,
models run in parallel and each model executes its CORDA stages serially. When
there are fewer models, models run sequentially and each stage distributes
`target direction × repeat` chunks across the full pool. A supplied inactive
`BPPARAM` is started once for the complete Layer-2 call and released at exit; an
already active caller pool remains active. Solvers are restricted to one thread
per worker to prevent oversubscription.

The persistent native solver removes repeated model construction and permits
basis reuse. A quantitative speedup claim requires a Human-GEM benchmark; the
synthetic checks establish correctness, not production-scale performance.

## Pipeline contracts

- Only retained HC core directions that can reach `corda_epsilon` enter the
  downstream `vmax`, penalty and score matrices.
- MC/NC/OT reconstruction does not enlarge the score matrix.
- Conditions and matching metacells share one model only within the same cell
  type and medium.
- The primary multiome path and RNA-only control reuse the same CORDA model file
  and checksum.
- CORDA caches are isolated under `model_cache/meta_module_gem/corda`.
- `fastcore_epsilon`, FASTCC counts, role-blocking counts and
  `max_support_reactions` are not reported as CORDA algorithm results.

## Audit output

Each model records the initial confidence, final RE status, inclusion stage,
reaction evidence, dependency tasks, Stage-2 NC-to-MC association pairs,
opposite-target directions blocked during each solve, solver backend, persistent
solver fallback counts, unpruned parent-model contract and HC-direction closure
checks. The full task table is stored once and a compact stage/status/backend
summary is stored separately.
