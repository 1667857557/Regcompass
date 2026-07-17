# Biological and mathematical correctness corrections

This note documents the canonical RegCompassR model contract used by v1.4.3. The objective is not to claim quantitative intracellular flux from RNA and ATAC. The primary output remains a medium- and GEM-dependent **minimum evidence-discordance penalty**.

## Implementation stages

The base implementations are refined by five explicitly ordered source files:

1. `workflow_stage_01_architecture.R` establishes the biological and inference contracts;
2. `workflow_stage_02_compatibility.R` applies compatibility refinements;
3. `workflow_stage_03_signed_projection.R` enforces signed-projection behavior;
4. `workflow_stage_04_result_contracts.R` enforces normalization, medium, penalty and result contracts; and
5. `workflow_stage_05_api_contracts.R` enforces the final public API and diagnostics contracts.

Their load order is declared in `DESCRIPTION` through `Collate`; it no longer
depends on opaque filenames containing increasing numbers of `z` characters.

## 1. Absolute evidence and relative state are separate

Metacell RNA input is `log1p(CPM)`. Canonical bounded support is computed from the supplied non-negative normalized signal:

\[
A_{gm}=\frac{x_{gm}}{x_{gm}+\kappa},
\qquad x_{gm}=\log(1+\operatorname{CPM}_{gm}).
\]

For ATAC and TF activity, the same zero-preserving transform is applied to the supplied non-negative normalized signal on that assay's own scale. The default is \(\kappa=1\), and assay-specific thresholds are configurable and recorded with the run.

\[
x_{gm}=0 \Rightarrow A_{gm}=0.
\]

A constant positive signal remains positive, whereas an all-zero gene remains zero. `A` is a bounded support score, not a probability, enzyme concentration or physical capacity. The previous within-gene robust-z sigmoid remains available only as a relative-state diagnostic and is `NA` for constant rows.

Reaction-wise Q95 scaling is retained as `C_within_reaction_relative`, but it is diagnostic only. `C_rel` is retained as a compatibility field and contains the same bounded absolute reaction evidence as `C_abs`; it is not a Q95-relative LP capacity.

## 2. Conservative GPR semantics

Canonical calls use explicit function defaults:

- required subunits: `AND = min`;
- alternative isoenzymes: `OR = max`;
- promiscuity correction: `none`.

These defaults avoid soft compensation of missing required subunits, isoenzyme-count inflation from summed OR rules, and annotation-degree-dependent penalties from promiscuity correction. Boltzmann AND, summed OR and promiscuity heuristics remain available only when callers request them explicitly for sensitivity analysis. The result no longer depends on the global `RegCompassR.strict_gpr_defaults` option.

A recursive Boolean parser supports nested AND/OR expressions and converts them to disjunctive normal form. Human-GEM import stops with reaction context on any malformed rule; parser errors are not silently converted to missing GPRs. Long-table GPR input must contain `and_group_id`.

## 3. Signed Pando regulatory support

Pando evidence is computed from signed TF-peak-gene activity:

\[
u_{em}=\operatorname{sign}(\hat\beta_e)T_{tm}A_{pm},
\]

weighted by \(|\hat\beta_e|\sqrt{R_e^2}\). Regulatory support is centered at 0.5:

\[
F_{gm}=\operatorname{clip}_{[0,1]}
\left(0.5+0.5\sum_e w_eu_{em}\right).
\]

Thus 0.5 is neutral, values below 0.5 represent active repression, and values above 0.5 represent active support. TF RNA support and peak accessibility are both required for a non-neutral edge contribution. Missing regulatory evidence is neutral rather than equivalent to missing enzyme expression.

The shared-TF projection retains regulator sets, direct regulator/target fields, signed projection weights and concordant/discordant/mixed regulatory relations. Direct TF-to-target edges and concordant shared-TF edges may define components; discordant or mixed shared-TF relations remain diagnostic and do not merge genes into one biological module.

## 4. Multiome penalty

The effective v1.4.3 expression penalty is the bounded inverse-support term:

\[
P^{expr}_{rm}=1-C^{abs}_{rm}.
\]

Observed zero expression, missing expression and no-GPR expression evidence receive the same maximum expression penalty of 1. Missing expression does not receive a second additive missingness penalty.

Regulatory support is converted to a repression-only modifier:

\[
F^{eff}_{rm}=\min\{1,\max(2F_{rm},\epsilon)\},
\]

\[
P^{conf}_{rm}=-\log(F^{eff}_{rm}).
\]

Consequently, neutral or positive Pando support (\(F\geq0.5\)) adds no penalty, while active repression (\(F<0.5\)) adds a non-negative penalty. The total biological penalty is

\[
P_{rm}=w_E P^{expr}_{rm}+w_F P^{conf}_{rm}
+w_G GPR_{missing}+P^{role}_{r}.
\]

The default expression component is bounded by 1, whereas the repression component can be larger when \(F\) approaches zero. Therefore `confidence` weight is a model-strength parameter, not a cosmetic display setting. Analyses should report sensitivity to the Pando confidence weight, including an RNA-only setting with `confidence = 0`.

FASTCORE membership is structural support and does not itself prove expression or regulation. Exchange, demand, sink and other structural roles may receive explicit role penalties; these overrides must be reported separately from biological evidence components.

## 5. Full-transcriptome RNA normalization is required

Metacell logCPM for GPR genes must use library sizes computed from the full transcriptome before filtering to GPR genes. The canonical workflow caches those full-library sizes during metacell construction and consumes them during GPR normalization.

If the cache is unavailable, internal callers must pass `library_size` explicitly. RegCompassR no longer silently substitutes the GPR-subset column sums, because that fallback can substantially inflate metabolic-gene CPM values.

## 6. Medium semantics

The default `compass_model_bounds` scenario preserves Human-GEM exchange directionality and caps each exchange bound at a shared COMPASS-style limit. It is a **shared model-defined environment**, not a measured physiological medium.

Concentration is not converted to uptake rate. Custom media should provide measured or justified reaction bounds and their provenance. Control and perturbation groups should share the same medium unless external measurements support condition-specific exchange constraints.

## 7. Inference unit

The canonical inference unit is `sample_celltype`. Metacell scoring remains available for within-sample heterogeneity but emits a warning because metacells from the same animal or patient are not independent biological replicates.

## 8. Score semantics

The LP penalty is the primary output. The display score is a stable within-target empirical rank:

\[
R_{rm}=1-\frac{\operatorname{rank}(P_{rm})-1}{n_r-1}.
\]

It is not a probability and is not intended for comparison across unrelated target reactions. Constant targets are marked non-informative and receive `NA` rather than an epsilon-amplified pseudo-score.

## Validation invariants

The tests require:

1. all-zero genes produce zero absolute evidence;
2. constant positive genes retain positive absolute evidence;
3. nested GPR rules preserve Boolean logic;
4. default required-subunit capacity equals the bottleneck;
5. default isoenzyme capacity is not inflated by isoenzyme count;
6. GPR defaults do not depend on a hidden global option;
7. missing regulatory evidence is neutral and active repression is penalized;
8. structural support cannot silently receive an unreported biological interpretation;
9. GPR-subset normalization fails without full-transcriptome library sizes;
10. constant penalties yield no relative rank;
11. model-bound media preserve blocked exchanges and model directionality;
12. discordant shared-TF relations do not merge components;
13. sample-by-cell-type is the default inference unit.

## Remaining scope

These corrections improve internal validity but do not make RNA+ATAC sufficient for absolute flux identification. Quantitative flux claims still require uptake and secretion measurements, proteomics or enzyme constraints, and ideally matched isotope tracing. Pando coefficient signs should be treated as model-derived regulatory evidence, not causal proof, unless supported by perturbation or independent regulatory validation.
