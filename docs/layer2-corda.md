# Original CORDA2 in Layer 2

`rc_regcompass_step_layer2()` keeps compact add-only FASTCORE as the default. The optional CORDA2 route is available with `model_mode = "meta_module_gem"` and follows the original MATLAB implementation in `schultzdre/Constraint-Based-Modeling/CORDA2.m`.

## Parameters

```r
layer2_args = list(
  target_direction = "both",
  solver = "highs",
  model_params = list(
    model_completion = "corda2",
    completion_time_limit = 3000,
    strict = TRUE,
    corda2_args = list(
      MCxNCthresh = 2,
      constraint = 1,
      constrainby = "val",
      om = 1e4,
      ci = 0.01
    ),
    corda_medium_confidence_threshold = 0.75,
    corda_negative_confidence_threshold = 0.10,
    corda_regulatory_weight = 0.20,
    corda_include_evidence_outside_modules = TRUE,
    corda_max_medium_confidence_reactions = Inf
  )
)
```

| Argument | Original default | Effect |
|---|---:|---|
| `MCxNCthresh` | `2` | Minimum number of remaining MC directional reactions that must depend on an NC direction before that NC direction is retained. |
| `constraint` | `1` | Absolute target flux under `constrainby = "val"`, or percentage under `"perc"`. |
| `constrainby` | `"val"` | Target constraint mode: absolute value or percentage of maximum flux. |
| `om` | `1e4` | High cost assigned to NC reactions and Step-3 OT reactions; Step-1 MC cost is `sqrt(om)`. |
| `ci` | `0.01` | Proportional cost increase applied once to each newly observed high-cost dependency during one target assessment. |

## Confidence mapping

RegCompass maps cell-type evidence to the four reaction groups required by original CORDA2:

| RegCompass reaction class | CORDA2 group |
|---|---|
| merged core reaction | HC |
| non-core merged meta-module reaction and selected high-evidence outside-module reaction | MC |
| finite low-evidence reaction | NC |
| remaining reaction | OT |

This mapping is an input adapter. The reconstruction state machine begins after HC, MC, NC and OT have been assigned.

## Native Layer 1 evidence aggregation

The canonical Layer 1 output stores RNA-only reaction capacity, multiome reaction capacity and reaction-level regulatory support as `reaction × metacell` matrices. CORDA2 reconstructs one shared structural model per `cell type × medium`, so these metacell values are converted to one evidence value per `reaction × cell type` before HC, MC, NC and OT classification.

For each reaction, cell type and condition, RegCompass first takes the median across metacells:

\[
E_{r,c,k}^{\mathrm{within}}
=
\operatorname{median}_{u \in (c,k)} E_{r,u}.
\]

It then takes the maximum of those condition medians within the same cell type:

\[
E_{r,c}^{\mathrm{union}}
=
\max_k E_{r,c,k}^{\mathrm{within}}.
\]

This fixed `condition_median_max` policy is applied independently to:

- `reaction_expression_rna_only`;
- `reaction_expression`;
- `reaction_regulatory_support_fraction`, when available.

The within-condition median limits the influence of individual metacell outliers. The cross-condition maximum preserves reactions strongly supported in any retained condition so that condition-specific pathways are not averaged away from the shared union GEM. No new user parameter is introduced, and all existing CORDA2 thresholds and original MATLAB controls remain unchanged. Legacy reaction-support data frames remain accepted and retain their previous cell-type aggregation behavior.

## Directional model

Only reactions satisfying `lb < 0 && ub >= 0` are split. Their reverse copies use the original suffix `_CORDA_rev_rxn`. The original forward column remains unchanged and the reverse copy uses the negated stoichiometric column. Both variables are non-negative.

When a directional target is tested, its opposite directional copy is closed. The target is first maximized and then fixed according to `constraint` and `constrainby`:

- `"val"`: `min(constraint, maximum target flux)`;
- `"perc"`: `0.01 * constraint * maximum target flux`.

## Reconstruction steps

### Step 1

Each HC direction is tested with the cost vector:

\[
c_i =
\begin{cases}
\sqrt{om}, & i \in MC,\\
om, & i \in NC,\\
10^{-3}, & \text{otherwise}.
\end{cases}
\]

Newly observed MC or NC dependencies have their cost multiplied by `1 + ci`. Solving continues until the HC-to-MC and HC-to-NC dependency sets stop changing. Used MC and NC directions are promoted to HC. Blocked HC directions are removed unless they were active while supporting another HC target.

### Step 2.1

Each remaining MC direction is tested with NC cost `om`. The MC-to-NC dependency matrix is updated until stable. MC directions unable to carry the required flux are removed.

### Step 2.2

NC directions required by at least `MCxNCthresh` MC directions are promoted. All other NC directions are blocked. Each remaining MC direction is then maximized with its opposite direction closed; infeasible MC directions are removed and recorded in `rescue`.

### Step 3

All reactions outside retained HC and OT are blocked. Each retained HC direction is tested with OT cost `om`; required OT directions are added to the final model.

## Final reaction bounds

Forward and reverse directional selections are merged back to the original reaction identifiers. RegCompass retains the medium-constrained parent bounds and disables directions not selected by CORDA2. This preserves the requested medium rather than reopening exchanges after reconstruction.

## Engineering flow

One CORDA2 reconstruction is created per `cell type × medium` combination. Target assessments inside one reconstruction remain serial because dependency matrices and promoted reaction sets are updated in source order. Independent reconstructions may run through `BPPARAM`. The solver is restricted to one thread per worker, and the completed structural model is reused by the primary and RNA-only scoring paths.

Each model stores the HC-to-MC, HC-to-NC and MC-to-NC matrices, rescue table, directional inclusion set, solver diagnostics, parameter values and cache checksum.
