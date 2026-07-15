# RegCompassR 1.4.0 (development)

## Global metacell workflow

- Replaced the staged integrated workflow with one upstream worker per retained `condition × sample × cell-type` stratum. Each worker now completes metacell construction, fragment aggregation, LinkPeaks/Layer 1, Pando, and meta-module inference before returning an artifact.
- Added a hard all-strata/all-sample barrier. Global recalibration and GEM construction are blocked when any retained stratum fails or any biological sample is absent.
- Recomputed reaction-capacity normalization after the barrier using one reaction-wise Q95 across all metacells.
- Unioned all strict-stratum reaction envelopes into one canonical `GLOBAL_UNION` meta-module and completed one shared GEM per medium scenario.
- Required condition-invariant medium constraints for shared-GEM comparison.
- Changed Layer 2 parallelization to one shared-model/medium × metacell task, with one metacell-specific penalty vector and all target directions evaluated after loading the model once.
- Explicitly releases the upstream worker pool before global processing and creates a fresh Layer 2 pool; cleanup also runs on errors.

## API cleanup

- Removed the obsolete staged runner `rc_run_regcompass_multiome_metacell()`.
- Removed the sample-specific `rc_build_meta_module_gem_cache()` implementation.
- Removed `rc_load_metacell_object_from_run()` and the versioned `R/regcompass_v13.R` file.
- Corrected parent-GEM failure semantics so an allowed reaction direction that is structurally infeasible is reported as `parent_blocked`, not `no_allowed_direction`.

# RegCompassR 1.3.0

## Architecture

- Restricted Layer 2 structural modeling to exactly two modes: `full_gem` and `meta_module_gem`.
- Removed the public target-k-hop, module-meso-GEM, and automatic fallback interfaces.
- Added `rc_run_regcompass()` as the integrated Layer 1, Pando meta-module, and Layer 2 entry point.
- Added one completed-model cache per sample, GRN module, and medium scenario.

## FASTCORE

- Added medium-specific parent-GEM feasibility validation and FASTCC consistency screening.
- Added the FASTCORE LP-7 and LP-10 sequence for add-only support completion.
- Implemented the original LP-10 scaling convention: core constraints and flux bounds are multiplied by `1e5`, while support is extracted using the original epsilon.
- Preserved the complete biological reaction envelope and penalized only candidate support reactions outside that envelope.
- Added signed reaction orientation for reverse-only and reversible core tasks without splitting reversible reactions into artificial forward/reverse copies.
- Demand, sink, and artificial-support reactions are disabled during structural completion; parent-blocked targets are reported instead of gap-filled.

## microCOMPASS

- Replaced the loose positive/negative flux split with signed reaction variables and absolute-value auxiliary variables.
- Non-zero positive and negative reaction bounds are now preserved exactly in the penalty LP.
- Meta-module scoring always uses the completed local GEM's own directional `vmax`; no alternative reference modes are exposed.
- Sample-specific meta-modules are evaluated only against Layer 2 units from the matching biological sample.

## Documentation and validation

- Replaced the v1.2 module-meso tutorial with the implemented v1.3 full-GEM/meta-module-GEM workflow.
- Added mathematical design documentation for FASTCC, FASTCORE, directional flux scoring, cache scope, and failure semantics.
- Added regression tests for LP-7/LP-10 support completion, paper scaling, reverse-only targets, forced non-zero bounds, parent-blocked reactions, reversible-cycle avoidance, cache reuse, and the two-mode API contract.

# RegCompassR 1.2.0

## Added

- Sample-specific Pando GRN inference on retained RNA+ATAC metacells.
- Metabolic target-gene selection from the intersection of single-cell RNA features and Human-GEM GPR genes.
- Export of complete and significant TF–peak–gene coefficient tables for every sample.
- Projection of significant Pando edges into sample-specific metabolic gene networks using shared-TF and direct metabolic-TF relationships.
- GRN-defined core reaction mapping and ordered meta-module expansion through core subsystems, shared KEGG/Reactome identifiers at subsystem level, and shared master-Rhea identifiers at subsystem level.
- Optional fixed-point expansion for sensitivity analysis.
- Human-GEM annotation preparation retaining subsystem, KEGG, Reactome, Rhea, and master-Rhea fields.

# RegCompassR 1.1.0

- Added strict condition × sample × cell type metacell filtering gates before and after metacell construction.
- Hardened Human-GEM archive download fallback and ZIP validation.
- Updated the formal metacell README example for the current API.

## Mathematical and biological audit

- Enforced exact penalty/reaction alignment and correct ranged-constraint handling.
- Added condition-aware full-GEM and meta-module parent caches.
- Reported blocked core reactions as `no_allowed_direction`.
- Marked FASTCORE as an add-only direction-preserving extension.
- Moved all differential inference to biological-sample aggregates.
- Activated documented GPR/role penalty parameters and made artificial drains expensive by default.
