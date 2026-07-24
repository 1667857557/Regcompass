# RegCompassR 1.8.1

- Added formal reaction annotation to Stage 6 and condition-statistics outputs: reaction names, stoichiometry-derived formulas with metabolite names and compartments, direction-specific substrates/products, subsystems, GPR rules, participating genes, and database cross-references.
- Added condition-by-cell-type evidence provenance that distinguishes active `RNA+ATAC` support from `RNA-only`, `GPR/no-observed-RNA`, and `structural/no-GPR` reactions. `RNA+ATAC` now requires the GPR-aggregated reaction capacity calculated from integrated evidence to differ from the otherwise identical RNA-only reaction capacity; gene-level ATAC modifiers and contribution genes are reported separately.
- Added `rc_build_reaction_annotations()` and `rc_attach_reaction_annotations()` for new and previously generated results.
- Added `rc_select_gene_reactions()` and `rc_plot_condition_gene_reactions()` for selecting scored reactions by metabolic genes and generating a ranked collection of significant, biologically annotated condition boxplots.
- Added `rc_test_condition_reactions()` for same-reaction, same-direction, same-medium comparisons between conditions within each cell type under the shared union-GEM. It reports Kruskal-Wallis omnibus tests, pairwise Wilcoxon tests, BH-adjusted P values, median score shifts, rank-biserial/common-language effects, and Cohen's d.
- Added `rc_plot_condition_reaction()` for multi-condition boxplots of a selected reaction target, with every metacell shown as a jittered point, Kruskal-Wallis omnibus annotation, and pairwise significance brackets based on raw or reaction-wide multiplicity-adjusted P values.
- Condition-reaction statistics explicitly distinguish within-dataset metacell significance from biological-replicate-level treatment inference and verify that target `vmax` is invariant across units before testing.
- Fixed `mouse_plasma` so it no longer inherits human HPLM concentrations or provenance. Healthy-mouse glucose (4.381 mM), lactate (3.088 mM), and glutamine (0.934 mM) define the only quantitative relative uptake caps; all other mouse components are availability-only.
- Separated the healthy-mouse quantitative reference from the broader murine plasma and tumor-interstitial-fluid availability evidence, and removed the unrelated Mouse-GEM reconstruction DOI from medium-composition provenance.
- Removed the redundant public `metacell_label_col` and stepwise `label_col` arguments. The canonical workflow now exposes its actual behavior directly: `celltype_col` is always passed to SuperCell2 before aggregation, while condition remains the only hard metacell stratum.
- Retained `label_col` only on the lower-level general-purpose `rc_make_supercell2_metacells()` builder, where it is a functional SuperCell2 option.
- Updated the README, all three tutorial levels, the workflow vignette, API index, and help pages to use the canonical interface only.
- Added a complete guide to the predefined extracellular media, including species restrictions, culture versus plasma backgrounds, glucose/lactate/glutamine sensitivity bounds, technical baselines, custom media, and the rule that medium constraints never expand original GEM directionality.

# RegCompassR 1.7.0

- Changed the canonical metacell scope to `condition × cell type`, deliberately pooling cells from all biological samples within each condition before SuperCell2 while retaining per-metacell biological-sample composition diagnostics.
- Changed Pando inference and GRN meta-module construction to the same condition-by-cell-type scope.
- Allows Pando installed from a locally downloaded source archive when GitHub remote metadata are unavailable. Such installations continue with an explicit warning and are marked as having an unverified repository origin; explicitly conflicting remote username or repository metadata still fail.
- Uses condition-specific Pando coefficients learned from RNA+ATAC to weight accessibility-only regulatory deviations at the metacell level; metacell TF RNA is not multiplied into the modifier, reducing direct duplicate RNA weighting.
- Clarifies that coefficients estimated from the same pooled dataset are fitted parameters rather than independent validation evidence; condition-pooled outputs remain descriptive unless external fitting or cross-fitting is supplied.
- Fixed the canonical GPR calculation to a normalized, monotone Boltzmann soft-min AND, additive isozyme OR, and no promiscuity weighting.
- Replaced the previous decomposed expression-plus-confidence objective with one COMPASS-like positive cost, `1 / (1 + log2(1 + E_multiome))`.
- Restricts fixed structural penalties to exchange, demand, sink, and artificial-support reactions. Transport and cofactor reactions with GPR evidence retain the integrated multiome reaction-expression cost.
- Builds biological meta-modules only from complete-GPR core reactions, core-reaction subsystems, and reactions sharing KEGG, Reactome, or master-Rhea identifiers. Metabolite-neighbour expansion is not used; local FASTCORE is the sole mechanism for adding reactions required for flux feasibility.
- Supports both shared union meta-module GEM and shared full-GEM scoring modes with the same Layer 1 evidence, medium, target-flux fraction, and ranking outputs.
- Allows one or more biological samples per condition. Sample counts are retained as provenance and do not block the descriptive pooled-metacell workflow.
- Allows one condition. Single-condition runs return within-condition reaction priorities; multi-condition runs additionally return all pairwise descriptive priority contrasts within each cell type.
- Added explicit `reaction_ranking` output containing reaction ID, direction, medium, median minimum penalty, support score, and within-condition priority rank.
- Deleted obsolete sample-level differential/statistics code and unused pseudobulk interfaces that were incompatible with the pooled-metacell inference semantics.
- Deleted the retired strict-stratum global workflow, Q95 calibration implementation, Pando reaction-confidence implementation, Layer 2 confidence alignment functions, confidence placeholders, `penalty_weights` API, and metabolite-neighbour expansion helper and controls.

# RegCompassR 1.6.0

- Added `fragment_files = FALSE` support so one-shot and integrated workflows can skip fragment aggregation and use object ATAC peak raw counts when matching fragment files are unavailable.
- Removed deprecated one-shot `humangem_version` handling; use `gem_version` with `species`.
- Clarified and tested that `rc_make_medium_scenarios()` can return preset and user-defined custom scenarios together, while preserving literature-derived concentration provenance and relative uptake sensitivity bounds.
- Updated tutorials, help pages, and public-API tests to document the canonical interfaces only.

# RegCompassR 1.4.2

- Fixed metacell RNA normalization so GPR-gene logCPM uses the full-transcriptome library size computed before filtering to metabolic genes.
- Replaced the expression term with a COMPASS-like inverse-support penalty, preventing missing/no-GPR evidence from receiving a lower penalty than observed zero expression.
- Added a shared `compass_model_bounds` medium that preserves GEM exchange directionality and caps exchange fluxes at a uniform limit of 1 by default.
- Applied structural penalties to exchange, demand, sink and artificial-support reactions independently of how their roles were annotated.
- Preserved the existing strict-stratum Pando workflow: peak-gene links remain inferred independently within each condition × sample × cell-type group.

# RegCompassR 1.4.1

- Replaced the canonical relative-z/Q95 LP capacity with zero-preserving absolute RNA evidence; Q95 is diagnostic only.
- Changed canonical integrated GPR defaults to hard-min AND, max OR and no promiscuity down-weighting.
- Added recursive nested Boolean GPR parsing and fail-fast Human-GEM import diagnostics.
- Reworked Pando evidence as signed TF–peak–gene regulatory support with TF expression and peak accessibility.
- Preserved regulator and sign metadata in shared-TF projections.
- Made missing/neutral regulatory evidence neutral in the LP penalty and prevented silent structural-support penalty overrides.
- Added explicit named medium backgrounds without retaining compatibility aliases for retired names.
- Changed the canonical inference unit to sample by cell type; metacell-level scoring is explicitly exploratory.
- Replaced the MAD-sigmoid display score with a stable within-target empirical penalty rank; raw penalty is the primary output.

# RegCompassR 1.4.0

- Focused the public API on the canonical workflow and its required setup helpers: `rc_prepare_human2_gem()`, `rc_make_medium_scenarios()` and `rc_run_regcompass()`.
- Kept tutorials concise while still showing adjustable setup steps for Human-GEM preparation and shared medium construction.
- Removed the adaptive metacell gamma API; the workflow now uses one fixed gamma and skips strata that do not produce enough metacells for downstream analysis.
- Removed standalone LinkPeaks, staged Layer 1, versioned Human-GEM and legacy reporting interfaces from the supported API surface.
