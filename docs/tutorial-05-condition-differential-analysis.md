# Tutorial Level 5: condition comparison from regulatory edge to directional reaction score

**Previous:** [Tutorial 4 — targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md) is optional and adds a second scoring pass for directly cross-referenced non-core reactions.

**This tutorial:** interprets the original Stage 5 scores by tracing one biological hypothesis through the complete RegCompass evidence chain:

```text
bootstrap-active TF–peak–target edge
→ condition-specific metabolic target gene
→ complete GPR branch
→ condition core reaction
→ shared medium-specific union GEM
→ RNA+ATAC reaction penalty
→ direction-specific LP support
→ descriptive condition contrast
```

**Inputs:** a compact `result` produced by Tutorial 1 or Stage 6 of Tutorial 2. Detailed Stage 1–5 objects are loaded only when a diagnostic requires them.

**Next:** use the exported evidence bundle for pathway interpretation, orthogonal validation, or a prespecified experimental follow-up. RegCompass does not convert metacell separation into biological-replicate inference.

---

## 1. Load the compact result and optional stage checkpoints

```r
result <- readRDS(
  "RegCompass_steps/06_results/regcompass_result.rds"
)

# Optional detailed objects. They are not duplicated inside `result`.
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

The compact result contains primary tables:

```r
result$table_manifest
names(result)
```

Detailed candidate networks, all coefficients, metacell matrices, full module membership, and Layer 1 matrices remain in the corresponding stage RDS and TSV files. This avoids storing the same high-dimensional information repeatedly.

Set the metadata and biological comparison once:

```r
condition_col <- "Group"
celltype_col <- "cell_type"

condition_a <- "Control"
condition_b <- "JQ1"
cell_type <- "stem-cell_like"
medium_id <- "high_glucose"
```

---

## 2. Confirm that the comparison is structurally valid

A condition comparison is interpretable only when all compared units use the same medium-specific structural model. Conditions may have different evidence penalties, but not different reaction sets, stoichiometry, or flux bounds.

```r
result$params$structural_comparability
result$microcompass$model_cache_summary
```

Inspect one row per medium:

```r
cache <- result$microcompass$model_cache_summary
cache[, intersect(c(
  "medium_scenario",
  "file",
  "file_checksum",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "build_strategy"
), colnames(cache)), drop = FALSE]
```

The comparison object also verifies that target `vmax` is invariant across metacells under the shared model. A failure here is a structural-model error, not a biological result.

The directional LP is:

\[
V_{r,d}^{\max}=\max_v d\,v_r
\quad\text{subject to}\quad
Sv=0,\;lb\le v\le ub,
\]

followed by the evidence-weighted near-optimal problem:

\[
\min_{v,z}\sum_j p_jz_j,
\quad z_j\ge |v_j|,
\quad d\,v_r\ge\omega V_{r,d}^{\max}.
\]

Because `S`, `lb`, `ub`, the medium, and `omega` are shared, differences between conditions arise from evidence-derived penalties rather than model topology.

---

## 3. Choose a reaction hypothesis without selecting by the same P value

For discovery, inspect the compact condition contrast and rank by effect magnitude or biological relevance. For confirmatory use, specify reactions before testing.

```r
contrast <- result$condition_contrast

candidate_contrasts <- contrast[
  contrast$cell_type == cell_type &
    contrast$medium_scenario == medium_id &
    contrast$condition_a == condition_a &
    contrast$condition_b == condition_b,
  , drop = FALSE
]

candidate_contrasts <- candidate_contrasts[
  order(-abs(candidate_contrasts$delta_support_b_minus_a)),
  , drop = FALSE
]

head(candidate_contrasts, 20)
```

Select a fixed reaction and direction:

```r
reaction_id <- "MAR04324"
direction <- "forward"
```

Retrieve the unique reaction annotation instead of repeating formula and GPR columns in every ranking row:

```r
reaction_annotation <- result$reaction_catalog[
  result$reaction_catalog$reaction_id == reaction_id,
  , drop = FALSE
]
reaction_annotation
```

---

## 4. Trace the condition-specific regulatory evidence

### 4.1 Bootstrap-active TF–peak–target edges

```r
edges <- result$active_regulatory_edges

edge_trace <- edges[
  edges[[condition_col]] %in% c(condition_a, condition_b) &
    edges[[celltype_col]] == cell_type,
  , drop = FALSE
]

edge_trace <- edge_trace[
  order(
    edge_trace[[condition_col]],
    -abs(edge_trace$stable_estimate)
  ),
  , drop = FALSE
]
head(edge_trace, 30)
```

For edge `e` and condition `c`, the fitted coefficient is

\[
\widehat\theta_{e,c}=\widehat\beta_e+\widehat\delta_{e,c},
\qquad \sum_c\widehat\delta_{e,c}=0.
\]

Bootstrap stability is

\[
\Pi_{e,c}=\frac{1}{B_s}\sum_{b\in\mathcal B_s}
I\!\left(|\widehat\theta^{(b)}_{e,c}|>\varepsilon\right),
\]

and conditional sign stability is

\[
\rho_{e,c}=\left|
\frac{\sum_b I_e^{(b)}\operatorname{sign}(\widehat\theta^{(b)}_{e,c})}
     {\sum_b I_e^{(b)}}
\right|.
\]

Layer 1 uses the reliability-weighted coefficient

\[
\widetilde\theta_{e,c}
=\widehat\theta_{e,c}\Pi_{e,c}\rho_{e,c}.
\]

Before interpreting an edge, check:

```r
summary(edge_trace$selection_frequency)
summary(edge_trace$sign_stability)
summary(edge_trace$bootstrap_success_fraction)
summary(edge_trace$cv_rsq)
```

Bootstrap measures stability under cell resampling. It does **not** create biological replicates or establish TF causality.

### 4.2 Condition-specific regulated metabolic genes

```r
targets <- result$condition_target_genes

target_trace <- targets[
  targets[[condition_col]] %in% c(condition_a, condition_b) &
    targets[[celltype_col]] == cell_type,
  , drop = FALSE
]
target_trace
```

Positive and negative active edges both establish that a target gene is regulated. The sign affects the ATAC modifier; gene membership itself is not restricted to activation.

### 4.3 Complete-GPR core reactions

```r
core <- result$core_reactions

core_trace <- core[
  core$reaction_id == reaction_id &
    core[[condition_col]] %in% c(condition_a, condition_b) &
    core[[celltype_col]] == cell_type,
  , drop = FALSE
]
core_trace
```

For reaction `r`, GPR branch `k`, and condition-specific target-gene set `G_c`:

\[
Core_{r,c}=1
\iff
\exists k:\;B_{r,k}\subseteq G_c.
\]

An incomplete enzyme complex cannot become a core merely because another isozyme branch is complete. The detailed branch-level audit remains available in:

```r
step3$condition_modules$core_gene_reaction[
  step3$condition_modules$core_gene_reaction$reaction_id == reaction_id,
  , drop = FALSE
]
```

A reaction may still be scored even when it is not a core in one condition, because Stage 5 uses the union of condition/cell-type biological modules and reuses one shared model. Core membership defines module construction provenance, not condition-specific model deletion.

---

## 5. Inspect RNA versus RNA+ATAC evidence

```r
evidence <- result$reaction_evidence

reaction_evidence <- evidence[
  evidence$reaction_id == reaction_id &
    evidence$condition %in% c(condition_a, condition_b) &
    evidence$cell_type == cell_type,
  , drop = FALSE
]
reaction_evidence
```

Interpret `evidence_class` as:

| Class | Meaning |
|---|---|
| `RNA+ATAC` | ATAC integration changes the GPR-aggregated reaction capacity. |
| `RNA-only` | RNA supports the reaction, but ATAC does not change reaction capacity above tolerance. |
| `GPR/no-observed-RNA` | The reaction has a GPR but no observed RNA-supported capacity. |
| `structural/no-GPR` | Structural reaction without a gene rule. |

The gene-level log-odds update is

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
     {1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

When no active condition edge exists, `R = 0`, so

\[
C^{MO}_{g,u}=C^{RNA}_{g,u}.
\]

For detailed metacell matrices:

```r
step4$gene_support_rna
step4$gene_regulatory_modifier
step4$gene_support_multiome
step4$reaction_expression
```

---

## 6. Read the within-condition ranking and descriptive contrast

```r
ranking <- result$reaction_ranking

one_ranking <- ranking[
  ranking$reaction_id == reaction_id &
    ranking$target_direction == direction &
    ranking$medium_scenario == medium_id &
    ranking$cell_type == cell_type &
    ranking$condition %in% c(condition_a, condition_b),
  , drop = FALSE
]
one_ranking
```

The normalized penalty is

\[
P^{norm}_{r,d,u}
=\frac{P_{r,d,u}}{\omega V^{\max}_{r,d}},
\]

and the reported support score is

\[
S_{r,d,u}=-\log\left(P^{norm}_{r,d,u}+\varepsilon\right).
\]

Therefore:

- lower normalized penalty means less evidence cost is required;
- higher support score means stronger support for the specified reaction direction;
- the score is a model-supported reaction potential, not a measured flux.

Inspect the precomputed descriptive contrast:

```r
one_contrast <- result$condition_contrast[
  result$condition_contrast$reaction_id == reaction_id &
    result$condition_contrast$target_direction == direction &
    result$condition_contrast$medium_scenario == medium_id &
    result$condition_contrast$cell_type == cell_type &
    result$condition_contrast$condition_a == condition_a &
    result$condition_contrast$condition_b == condition_b,
  , drop = FALSE
]
one_contrast
```

`delta_support_b_minus_a > 0` means greater support in `condition_b`.

---

## 7. Run metacell-level condition statistics

```r
statistics <- rc_test_condition_reactions(
  result,
  condition_col = condition_col,
  celltype_col = celltype_col,
  conditions = c(condition_a, condition_b),
  comparisons = list(c(condition_a, condition_b)),
  reaction_ids = reaction_id,
  target_directions = direction,
  medium_scenarios = medium_id,
  cell_types = cell_type,
  min_units = 5L,
  include_omnibus = FALSE,
  p_adjust_method = "BH",
  p_adjust_scope = "celltype_contrast_medium",
  include_scores = TRUE
)

statistics$pairwise
```

Use the outputs in this order:

1. `delta_median_score_b_minus_a` — direction and magnitude of the median shift;
2. `rank_biserial_b_minus_a` — ordinal effect size;
3. `common_language_b_greater_a` — probability that a randomly selected `b` metacell has greater support than an `a` metacell;
4. `p_adj` — within-dataset metacell separation after the declared multiple-testing family.

The metacells are derived from the same biological dataset. Their P values are **not** biological-replicate treatment inference. Without independent samples, report them as descriptive within-dataset association statistics.

Do not first select reactions by the smallest `p_adj` and then describe the same P value as confirmatory evidence.

---

## 8. Report reversible reactions without double counting

Forward and reverse are separate counterfactual LP targets. They are not positive and negative measurements of one net flux.

```r
direction_report <- rc_report_condition_directions(
  result,
  condition_col = condition_col,
  celltype_col = celltype_col,
  conditions = c(condition_a, condition_b),
  comparisons = list(c(condition_a, condition_b)),
  reaction_ids = reaction_id,
  medium_scenarios = medium_id,
  cell_types = cell_type,
  min_units = 5L,
  include_unit_metrics = TRUE,
  source_label = "original_stage5_shared_union_gem",
  outdir = "RegCompass_result/07_direction_report"
)
```

Primary COMPASS-style outputs remain direction specific:

```r
direction_report$directional_pairwise
```

Non-additive summaries are:

\[
S^{any}=\max(S^{forward},S^{reverse}),
\]

\[
B^{direction}=S^{forward}-S^{reverse}.
\]

```r
direction_report$reaction_pairwise
direction_report$direction_diagnostics[, c(
  "reaction_id", "condition", "direction_pair_status",
  "max_abs_forward_reverse_difference",
  "directionally_indistinguishable", "preferred_direction"
)]
```

`directional_balance` is support asymmetry, not net flux. When forward and reverse are indistinguishable, report that limitation instead of choosing a direction from numerical noise.

---

## 9. Plot the fixed reaction direction

```r
rc_plot_condition_reaction(
  result,
  reaction_id = reaction_id,
  target_direction = direction,
  medium_scenario = medium_id,
  cell_type = cell_type,
  condition_col = condition_col
)
```

Plot all metacells, medians, and effect direction. Do not label a visual separation as biological-replicate significance.

---

## 10. Quality-control checklist

Before biological interpretation, require all relevant checks to pass:

```r
qc <- list(
  model_shared = result$params$structural_comparability,
  grn_metacell_coverage = result$grn_metacell_group_coverage[
    result$grn_metacell_group_coverage[[celltype_col]] == cell_type,
    , drop = FALSE
  ],
  active_edges = edge_trace,
  target_genes = target_trace,
  core_reaction = core_trace,
  reaction_evidence = reaction_evidence,
  ranking = one_ranking,
  contrast = one_contrast,
  statistics = statistics$pairwise,
  direction_diagnostics = direction_report$direction_diagnostics
)
```

Minimum questions:

- Was the same union GEM reused across conditions?
- Was the target scored in the same direction and medium?
- Were enough metacells available in both conditions?
- Did the target model pass CV and bootstrap-completion thresholds?
- Is the regulatory sign stable?
- Does a complete GPR branch support core membership?
- Does ATAC actually change reaction-level capacity, or is the reaction RNA-only?
- Are forward and reverse distinguishable?
- Is the claim descriptive, mechanistic, or independently validated?

---

## 11. Export one compact evidence bundle

```r
outdir <- "RegCompass_result/08_condition_evidence_bundle"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(edge_trace,
          file.path(outdir, "active_regulatory_edges.csv"), row.names = FALSE)
write.csv(target_trace,
          file.path(outdir, "condition_target_genes.csv"), row.names = FALSE)
write.csv(core_trace,
          file.path(outdir, "complete_gpr_core_reaction.csv"), row.names = FALSE)
write.csv(reaction_annotation,
          file.path(outdir, "reaction_annotation.csv"), row.names = FALSE)
write.csv(reaction_evidence,
          file.path(outdir, "reaction_evidence.csv"), row.names = FALSE)
write.csv(one_ranking,
          file.path(outdir, "condition_ranking.csv"), row.names = FALSE)
write.csv(one_contrast,
          file.path(outdir, "descriptive_contrast.csv"), row.names = FALSE)
write.csv(statistics$pairwise,
          file.path(outdir, "metacell_statistics.csv"), row.names = FALSE)
write.csv(direction_report$direction_diagnostics,
          file.path(outdir, "direction_diagnostics.csv"), row.names = FALSE)
```

This bundle preserves the complete interpretive chain without exporting all structural candidates, every coefficient, or every intermediate matrix.

For non-core reactions discovered through Tutorial 4, run the same sequence on the target-union result and set a distinct `source_label`. Do not combine original and second-pass testing families unless the multiple-testing scope is explicitly redefined.
