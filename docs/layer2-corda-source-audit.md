# Python CORDA2 source-to-R audit

This audit fixes the reference implementation to `resendislab/corda` commit
`c02e06d50606bf93f23d8f2e6d6ade0e996ca70e`. The parity claim is restricted
to `met_prod = NULL`, because RegCompass Layer 2 does not create the Python
class's temporary mock metabolite-production reactions.

The audit distinguishes three layers:

1. **CORDA2 reconstruction semantics** — reproduced from the pinned Python
   class without mathematical correction;
2. **RegCompass confidence construction** — converts cell-type multiome
   evidence into the five confidence values expected by CORDA2;
3. **RegCompass post-reconstruction scoring** — filters retained core reaction
   directions for downstream `vmax`, penalty and score calculation. This does
   not mutate the CORDA2 confidence updates.

## Source correspondence

| Python source operation | Mathematical or state effect | RegCompass implementation | Parity rule |
|---|---|---|---|
| module constants `UPPER` and `CI` | `UPPER = 10^6`; redundancy multiplier `CI = 1.01` | `.rc_layer2_corda_options()`, `.rc_corda_split_model()`, `.rc_corda2_associated()` | Fixed, not configurable |
| constructor `self.tflux = 1` | target lower bound and Stage-2 comparison use `1` | `.rc_layer2_corda_options()`, `.rc_corda_target_bounds()`, `.rc_corda2_minimize_medium_targets()` | Fixed, not configurable |
| constructor stores solver feasibility tolerance | all direction-open, activity and association tests use the active solver tolerance | `.rc_corda2_solver_feasibility_tolerance()`, `split$tolerance` | Solver-specific default is recorded in every model |
| constructor normalizes reaction bounds | lower bounds below `-tol` become `-10^6`; upper bounds above `tol` become `10^6` | `.rc_corda_split_model()` | Assignment order and threshold inequalities are preserved |
| COBRA reaction variables | every reaction has forward and reverse nonnegative variables, including a closed copy | `.rc_corda_split_model()` | Variable order is reaction order, forward then reverse |
| confidence expansion | one reaction confidence is copied to both direction variables | `.rc_corda2_directional_confidence()` | Five values only: `-1,0,1,2,3` |
| `associated()` penalty construction | the **forward variable confidence** chooses one cost, then that cost is assigned to both directions | `.rc_corda2_penalties()` | Direction-specific cost correction is deliberately not applied |
| `associated()` target bounds | if `ub < tol`, target becomes impossible; otherwise target `lb` is set first to `max(1,lb)`, then `ub` to `10^6` | `.rc_corda_target_bounds()`, `.rc_corda2_associated()` | Opposite reversible direction remains open; transient bound-order error is reproduced |
| dependency LP | minimize `c^T v` subject to `S_split v = 0` and current bounds | `.rc_corda_engine_solve()` | Full objective and bounds are reset before every solve |
| association test | include non-target variables with `v_i > tol` and confidence in `{-1,1,2}` | `.rc_corda2_associated()` | Uses solver primal values and solver tolerance |
| redundancy loop | newly observed penalized variables have cost multiplied by `1.01`; stop on no new variable or after `n` iterations | `.rc_corda2_associated()` | Iterations for one target are serial |
| target processing order | target list is processed serially on one mutable CORDA object and one solver state | `.rc_corda2_associated()` | No target-level parallelism inside a model |
| build iteration 1 | assess all confidence-3 directions with medium penalties enabled; promote associated directions to 3 | `.rc_corda_build_three_stage()` Stage 1 | Promotion occurs after the complete `associated()` call |
| build iteration 2a | assess remaining confidence-1/2 directions without medium penalties; count repeated absent associations; promote counts `>= support` | `.rc_corda_build_three_stage()` Stage 2 association/support | Counts retain duplicates across different signed targets |
| absent blocking | set each remaining confidence--1 variable upper bound to `max(0,lb)` | `.rc_corda_build_three_stage()` before Stage 2b | Direction variables are blocked individually |
| build iteration 2b | under a minimization objective, set coefficient `+1` for each remaining confidence-1/2 variable; promote only if the minimum objective is `> 1` | `.rc_corda2_minimize_medium_targets()` | The source's positive-coefficient minimization is preserved, not replaced by maximum flux |
| build iteration 3 preparation | block remaining confidence-1/2 variables with `ub = 0`; convert confidence 0 to -1 | `.rc_corda_build_three_stage()` Stage 3 preparation | Python setter failure for positive lower bounds is reproduced |
| build iteration 3 association | assess all confidence-3 targets once, absent penalties only, no redundancy; promote associated variables | `.rc_corda_build_three_stage()` Stage 3 | Serial target order remains unchanged |
| impossible targets | impossible variables are changed to -1 and retained in the final audit list | `.rc_corda2_associated()`, reconstruction output | Reconstruction completes; `strict` does not override source behavior |
| final inclusion | reaction is retained when `max(conf_forward, conf_reverse) == 3` | `.rc_corda2_reduce_confidence()` and `included` | Original medium-constrained bounds are restored in the final GEM |

## Mathematical programs

### Dependency assessment

For a signed target variable `x`, let `c` be the current penalty vector. The
source solves

\[
\begin{aligned}
\min_v\quad & c^T v\\
\text{s.t.}\quad &S_{split}v=0,\\
&l\le v\le u,\\
&v_x\ge 1.
\end{aligned}
\]

The source does **not** add `v_opposite = 0`. Therefore, an otherwise isolated
reversible target may satisfy the constraint through equal forward and reverse
flux. RegCompass preserves that behavior for exact parity.

For each reaction `r`, both directional coefficients are selected from the
forward confidence:

\[
c_{r,+}=c_{r,-}=\begin{cases}
1,&q_{r,+}\in\{1,2\}\text{ and medium penalties are enabled},\\
\text{penalty\_factor},&q_{r,+}=-1,\\
0,&\text{otherwise}.
\end{cases}
\]

After each optimum, a direction is associated when its primal flux is greater
than the active solver feasibility tolerance and its confidence is `-1`, `1`
or `2`. Newly observed penalized directions receive

\[
c_i\leftarrow 1.01\,c_i.
\]

### Stage-2 independent medium test

The pinned source does not maximize a medium direction. It sequentially solves

\[
\min_v\ v_x
\]

with the current bounds and promotes `x` only when the optimum is strictly
greater than `1`. This behavior is reproduced even though a maximum-flux test
would be a different and often more intuitive algorithm.

## Ordering and parallelism

Exactness requires preserving mutable state order within one CORDA2 instance:

- reaction order determines forward/reverse variable insertion order;
- targets are processed serially;
- redundancy iterations are serial;
- Stage-2 medium variables are processed serially;
- confidence mutations from an impossible target are visible to later build
  phases exactly as in Python.

RegCompass may parallelize independent `cell type × medium` instances, because
they do not share confidence vectors, bounds or solver state. It does not
parallelize targets inside one instance.

## Executable oracle

The workflow `CORDA2 exact-source checks` installs the pinned Python commit and
runs `tests/corda2-python-reference-check.py`. The generated oracle is consumed
by `tests/corda-synthetic-check.R`. The paired tests cover:

- fixed constants and constructor defaults;
- solver feasibility tolerance;
- two redundant dependency paths;
- the reversible forward/reverse self-cycle;
- absent-support counting across signed targets;
- forced versus free behavior in the positive-minimization Stage 2b;
- Stage-3 unknown-reaction completion;
- persistent HiGHS versus one-shot fallback on nondegenerate networks;
- installed-package execution with a live SnowParam pool while inner CORDA2
  execution remains serial.

Different LP solvers may choose different members of a degenerate optimal face.
The parity tests therefore use nondegenerate networks for exact output
comparison and separately test the source-level ordering, objective, bounds and
state-transition contracts.
