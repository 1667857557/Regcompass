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
