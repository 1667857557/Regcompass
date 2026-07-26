# Direction-aware condition reporting

`rc_report_condition_directions()` is the final reporting layer for reversible
reaction targets. It follows the useful part of the COMPASS reporting contract:
forward and reverse LP targets remain separate primary results. It does not
inherit the ambiguous practice of treating both rows as two independent
biological discoveries when their unit-level scores are numerically identical.

## Why a separate report is required

For a reversible GEM reaction, RegCompass solves two independent
counterfactual optimisation problems:

```text
forward:  v[r]  >= omega * Vmax_forward[r]
reverse: -v[r]  >= omega * Vmax_reverse[r]
```

Each problem minimises the same expression-linked absolute-flux cost. Forward
and reverse therefore can be both available and can produce identical support
scores when the bounds, network alternatives, and evidence costs are symmetric.
That does not mean both directions carry flux in one solution. It means the
available evidence cannot distinguish the preferred direction.

The support score tested between conditions remains:

```text
support = -log(penalty / (omega * vmax) + eps)
```

Higher values indicate stronger support for the specified LP target. They are
not measured fluxes.

## Run the report

```r
direction_report <- rc_report_condition_directions(
  result,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  cell_types = c("epithelial_like", "stem-cell_like"),
  comparisons = list(
    c("control_24hr", "JQ1_24hr"),
    c("control_24hr", "MS177_24hr"),
    c("JQ1_24hr", "MS177_24hr")
  ),
  medium_scenarios = "high_glucose",
  min_units = 5L,
  direction_tolerance = 1e-10,
  source_label = "original_layer2_core",
  outdir = "RegCompass_result/07_direction_report"
)
```

## Primary COMPASS-style directional results

```r
direction_report$directional_pairwise
direction_report$directional_omnibus
```

These tables retain one row per:

```text
reaction x direction x medium x cell type x contrast
```

Use `tested_formula` to identify the actual biochemical direction. A forward
and a reverse result are separate LP objectives, not two simultaneous fluxes.

## Non-additive reaction-level summaries

```r
direction_report$reaction_pairwise
direction_report$reaction_omnibus
```

The `report_metric` column has two values.

### `any_direction_support`

```text
max(forward_support, reverse_support)
```

If only one direction is available, that direction is retained. This is a
non-additive summary of the best-supported available direction. It avoids
double counting when forward and reverse are identical.

Use it when the question is:

> Does this reversible reaction have stronger condition-associated support in
> at least one allowed direction?

It is not a total flux or a sum of directions.

### `directional_balance`

```text
forward_support - reverse_support
```

Positive values favour the model-defined forward target; negative values favour
the reverse target. In pairwise output, `direction_shift_b_minus_a` reports
`toward_forward`, `toward_reverse`, or `no_shift`.

This is directional support asymmetry. It is not net flux because the two scores
come from different optimisation problems.

The two metric families receive separate multiplicity adjustment within the
selected `p_adjust_scope`.

## Direction diagnostics

```r
direction_report$direction_diagnostics
```

One row is returned per:

```text
reaction x medium x cell type x condition
```

Important fields are:

- `direction_pair_status`:
  - `forward_only`;
  - `reverse_only`;
  - `bidirectional_distinguishable`;
  - `bidirectional_indistinguishable`;
  - `bidirectional_no_paired_scores`;
- `max_abs_forward_reverse_difference`;
- `directionally_indistinguishable`;
- `preferred_direction`;
- `median_any_direction_support`;
- `median_directional_balance`.

A `bidirectional_indistinguishable` result means that forward and reverse have
the same finite/missing pattern and differ by no more than
`direction_tolerance` across the selected units. Report it once as reversible
reaction support, not twice as independent forward and reverse discoveries.

## Original cores and target-union reactions

Target-union does not rescore reactions that were already Layer 2 cores. Build
one report from the original Layer 2 object and one from the second-pass object:

```r
core_report <- rc_report_condition_directions(
  result,
  reaction_ids = all_original_core_ids,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  medium_scenarios = "high_glucose",
  source_label = "original_layer2_core"
)

expanded_report <- rc_report_condition_directions(
  targeted_for_stats,
  reaction_ids = expanded_noncore_ids,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  medium_scenarios = "high_glucose",
  source_label = "target_union_second_pass_noncore"
)
```

Combine matching tables while preserving `source_label`:

```r
all_directional_pairwise <- rbind(
  core_report$directional_pairwise,
  expanded_report$directional_pairwise
)

all_reaction_pairwise <- rbind(
  core_report$reaction_pairwise,
  expanded_report$reaction_pairwise
)

all_direction_diagnostics <- rbind(
  core_report$direction_diagnostics,
  expanded_report$direction_diagnostics
)
```

Because the two reports were tested separately, recalculate BH correction over
the combined testing family before using a joint FDR threshold. For directional
results, include `target_direction` as part of the feature identity. For
reaction-level results, correct `any_direction_support` and
`directional_balance` as separate metric families.

## Interpretation boundary

The report supports statements such as:

> MS177 metacells show stronger support for at least one allowed direction of
> the carnitine-shuttle reaction, while forward and reverse remain numerically
> indistinguishable.

It does not support statements such as:

> Forward and reverse flux both increased.

or:

> The difference between forward and reverse support is the net flux.

Net flux direction requires additional objective, thermodynamic, metabolomic,
exchange-flux, isotope-tracing, or flux-sampling information.
