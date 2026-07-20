# RegCompassR 1.7.0

- Changed the canonical metacell scope to `condition × cell type`, deliberately pooling cells from all biological samples within each condition before SuperCell2 while retaining per-metacell biological-sample composition diagnostics.
- Changed Pando inference and GRN meta-module construction to the same condition-by-cell-type scope.
- Uses condition-specific Pando coefficients learned from RNA+ATAC to weight accessibility-only regulatory deviations at the metacell level; metacell TF RNA is not multiplied into the modifier, reducing direct duplicate RNA weighting.
- Clarifies that coefficients estimated from the same pooled dataset are fitted parameters rather than independent validation evidence; condition-pooled outputs remain descriptive unless external fitting or cross-fitting is supplied.
- Fixed the canonical GPR calculation to a normalized, monotone Boltzmann soft-min AND, additive isozyme OR, and no promiscuity weighting.
- Replaced the previous decomposed expression-plus-confidence objective with one COMPASS-like positive cost, `1 / (1 + log2(1 + E_multiome))`.
- Builds biological meta-modules only from complete-GPR core reactions, core-reaction subsystems, and reactions sharing KEGG, Reactome, or master-Rhea identifiers. Metabolite-neighbour expansion is not used; local FASTCORE is the sole mechanism for adding reactions required for flux feasibility.
- Preserved one shared union-GEM, one shared medium, common bounds and directional two-step COMPASS-like LP scoring across all conditions.
- Added direct descriptive condition summaries and two-condition reaction-support contrasts within each cell type.
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
