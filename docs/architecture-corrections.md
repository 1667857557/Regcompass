# Biological and mathematical correctness corrections

This note documents the canonical RegCompassR corrections introduced in v1.4.1. The objective is not to claim quantitative intracellular flux from RNA and ATAC. The canonical output remains a medium- and GEM-dependent **minimum evidence-discordance penalty**.

## 1. Absolute evidence and relative state are separate

Metacell RNA input is log1p(CPM). Canonical bounded support is now computed from the supplied non-negative normalized signal:

\[
A_{gm}=\frac{x_{gm}}{x_{gm}+\kappa},
\qquad x_{gm}=\log(1+\operatorname{CPM}_{gm})
\]

for RNA, while ATAC uses its non-negative normalized accessibility signal on its own input scale. The default is \(\kappa=1\), and the threshold is configurable and recorded with the run. This transform preserves the absolute zero without applying an RNA-specific inverse transform to ATAC:

\[
x_{gm}=0 \Rightarrow A_{gm}=0.
\]

A constant positive signal remains positive, whereas an all-zero gene remains zero. `A` is a bounded support score, not a probability, enzyme concentration or physical capacity. The previous within-gene robust-z sigmoid remains available only as a relative-state diagnostic and is `NA` for constant rows.

Reaction-wise Q95 scaling is retained as `C_within_reaction_relative`, but it is diagnostic only. `C_rel` is retained as a compatibility field and now contains bounded absolute reaction evidence used by the LP.

## 2. Conservative GPR semantics

Canonical integrated runs use:

- required subunits: `AND = min`;
- alternative isoenzymes: `OR = max`;
- promiscuity correction: `none`.

These defaults avoid soft compensation of missing required subunits, isoenzyme-count inflation from `sum/sqrt(K)`, and annotation-degree-dependent penalties from `1/sqrt(K)`. Boltzmann AND, summed OR and promiscuity heuristics remain available to direct callers for sensitivity analysis.

A recursive Boolean parser now supports nested AND/OR expressions and converts them to disjunctive normal form. Human-GEM import stops with reaction context on any malformed rule; parser errors are no longer silently converted to missing GPRs. Long-table GPR input must contain `and_group_id`.

## 3. Signed Pando regulatory support

Pando evidence is now computed from signed TF–peak–gene activity:

\[
u_{em}=\operatorname{sign}(\hat\beta_e)T_{tm}A_{pm},
\]

weighted by \(|\hat\beta_e|\sqrt{R_e^2}\). Regulatory support is centered at 0.5:

\[
F_{gm}=\operatorname{clip}_{[0,1]}\left(0.5+0.5\sum_e w_eu_{em}\right).
\]

Thus 0.5 is neutral, values below 0.5 represent active repression, and values above 0.5 represent active support. TF RNA support and peak accessibility are both required for a non-neutral edge contribution. Missing regulatory evidence is neutral rather than equivalent to missing enzyme expression.

The shared-TF projection retains regulator sets, direct regulator/target fields, signed projection weights and concordant/discordant/mixed regulatory relations. Direct TF-to-target edges and concordant shared-TF edges may define components; discordant or mixed shared-TF relations remain diagnostic and do not merge genes into one biological module.

## 4. Multiome penalty

RNA remains the primary enzyme-availability evidence. Signed Pando support is a conservative modulation: missing or neutral regulation adds no penalty, while active repression adds a non-negative term.

\[
P_{rm}=w_E[-\log(C^{abs}_{rm})]
+w_F[-\log\{\min(1,\max(2F_{rm},\epsilon))\}]
+w_MI(\text{missing RNA})+w_GGPR_{missing}.
\]

FASTCORE membership is structural support and cannot silently overwrite this biological penalty. The legacy support-override helper now requires an explicit opt-in.

## 5. Medium semantics

The former generic `blood_like` baseline opened every exchange reaction with the same arbitrary uptake bound. It is now explicitly named `permissive_all_exchange` and labelled a technical upper-bound scenario. Legacy names emit warnings and are marked as non-curated sensitivity aliases.

Concentration is not converted to uptake rate. Custom media should provide measured or justified reaction bounds and their provenance. In biological comparisons, Control and perturbation groups should share the same medium unless external measurements support condition-specific exchange constraints.

## 6. Inference unit

The canonical inference unit is now `sample_celltype`. Metacell scoring remains available for within-sample heterogeneity but emits a warning because metacells from the same animal are not independent biological replicates.

## 7. Score semantics

The LP penalty is the primary output. The previous MAD-sigmoid display score is replaced by a stable within-target empirical rank:

\[
R_{rm}=1-\frac{\operatorname{rank}(P_{rm})-1}{n_r-1}.
\]

It is explicitly **not a probability**. Constant targets are marked non-informative and receive `NA` rather than an epsilon-amplified pseudo-score.

## Validation invariants

The added tests require:

1. all-zero genes produce zero absolute evidence;
2. constant positive genes retain positive absolute evidence;
3. nested GPR rules preserve Boolean logic;
4. canonical required-subunit capacity equals the bottleneck;
5. canonical isoenzyme capacity is not inflated by isoenzyme count;
6. missing regulatory evidence is neutral and active repression is penalized;
7. structural support cannot silently receive a low biological penalty;
8. constant penalties yield no relative rank;
9. permissive media are explicitly identified as technical assumptions;
10. shared-TF projection retains sign and regulator metadata;
11. discordant shared-TF relations do not merge components;
12. sample-by-cell-type is the default inference unit.

## Remaining scope

These corrections improve internal validity but do not make RNA+ATAC sufficient for absolute flux identification. Quantitative flux claims still require uptake and secretion measurements, proteomics or enzyme constraints, and ideally matched isotope tracing.
