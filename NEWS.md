# RegCompassR 1.4.1

- Replaced the canonical relative-z/Q95 LP capacity with zero-preserving absolute RNA evidence; Q95 is diagnostic only.
- Changed canonical integrated GPR defaults to hard-min AND, max OR and no promiscuity down-weighting.
- Added recursive nested Boolean GPR parsing and fail-fast Human-GEM import diagnostics.
- Reworked Pando evidence as signed TF–peak–gene regulatory support with TF expression and peak accessibility.
- Preserved regulator and sign metadata in shared-TF projections.
- Made missing/neutral regulatory evidence neutral in the LP penalty and prevented silent structural-support penalty overrides.
- Renamed the generic all-exchange medium to `permissive_all_exchange`; legacy biological-sounding aliases now warn.
- Changed the canonical inference unit to sample by cell type; metacell-level scoring is explicitly exploratory.
- Replaced the MAD-sigmoid display score with a stable within-target empirical penalty rank; raw penalty is the primary output.

# RegCompassR 1.4.0

- Focused the public API on the canonical workflow and its required setup helpers: `rc_prepare_human2_gem()`, `rc_make_medium_scenarios()` and `rc_run_regcompass()`.
- Kept tutorials concise while still showing adjustable setup steps for Human-GEM preparation and shared medium construction.
- Removed the adaptive metacell gamma API; the workflow now uses one fixed gamma and skips strata that do not produce enough metacells for downstream analysis.
- Removed standalone LinkPeaks, staged Layer 1, versioned Human-GEM and legacy reporting interfaces from the supported API surface.
