# Condition-associated reaction statistics

`rc_test_condition_reactions()` compares the same reaction target across
conditions within the same cell type after Layer 2 scoring.

## Why the comparison is valid

RegCompass builds one deduplicated global union-GEM from the locally completed
meta-modules identified across all condition-by-cell-type strata. Every metacell
is scored with the same stoichiometric matrix, bounds, medium, reaction
direction, and target-flux fraction. The only unit-specific term in the Layer 2
objective is the multiome-derived reaction penalty.

For reaction target `r` and unit `j`, the tested support score is

```text
S[r,j] = -log(P[r,j] / (omega * Vmax[r]) + eps)
```

where `P[r,j]` is the minimum evidence-discordance penalty required to sustain
`omega * Vmax[r]`. Larger scores indicate stronger support. Before testing, the
function verifies that `Vmax[r]` is invariant across units under the shared GEM.

## Basic use

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  cell_types = c("epithelial_like", "stem-cell_like"),
  min_units = 5,
  p_adjust_method = "BH",
  p_adjust_scope = "celltype_contrast_medium",
  outdir = "RegCompass_result/07_condition_statistics"
)
```

The function accepts either the complete RegCompass result or the Layer 2
`microcompass`/`step5` object.

## Outputs

`condition_stats$omnibus` contains one Kruskal-Wallis test per
`cell type × reaction × direction × medium` when at least three conditions are
available.

`condition_stats$pairwise` contains all requested pairwise Wilcoxon tests and:

- median and mean support scores in both conditions;
- `delta_median_score_b_minus_a`;
- Cohen's d;
- rank-biserial correlation;
- common-language probability that a condition-B unit exceeds condition A;
- raw and adjusted P values;
- the condition with higher predicted reaction support;
- explicit analysis-unit and inference-level labels.

Positive `delta_median_score_b_minus_a`, Cohen's d, or rank-biserial correlation
means stronger support in `condition_b`.

## Selected comparisons and reactions

```r
condition_stats <- rc_test_condition_reactions(
  step5,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  cell_types = "epithelial_like",
  comparisons = list(
    c("control_24hr", "JQ1_24hr"),
    c("control_24hr", "MS177_24hr"),
    c("JQ1_24hr", "MS177_24hr")
  ),
  reaction_ids = c("MAR06231", "MAR06241"),
  target_directions = c("forward", "reverse"),
  medium_scenarios = "high_glucose"
)
```

## Candidate filtering

P values should be interpreted together with effect sizes. A practical
exploratory filter is:

```r
hits <- subset(
  condition_stats$pairwise,
  p_adj < 0.05 &
    abs(rank_biserial_b_minus_a) >= 0.30 &
    abs(delta_median_score_b_minus_a) >= 0.10
)
```

Forward and reverse targets remain separate because they are separate LP
objectives. A significant score difference indicates differential support for a
reaction direction, not direct observation of net flux.

## Inference boundary

For `unit = "metacell"`, the P values quantify within-dataset separation of
condition-associated metacell score distributions. They do not make metacells
independent biological replicates. Results therefore carry:

```text
inference_level = metacell_within_dataset
descriptive_only = TRUE
biological_replicate_inference = FALSE
```

With independent biological samples, formal treatment-level inference should
use sample-by-cell-type units or a sample-aware hierarchical analysis.
