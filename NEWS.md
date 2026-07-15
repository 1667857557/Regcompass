# RegCompassR 1.4.0 (development)

- Removed the duplicate Signac LinkPeaks pass from the integrated workflow; Pando is now the sole peak-gene model, and significant Pando regions plus metacell accessibility define ATAC confidence.

## Global metacell workflow

- Replaced the staged integrated workflow with one upstream worker per retained `condition × sample × cell-type` stratum. Each worker now completes metacell construction, fragment aggregation, Pando peak-gene/GRN inference, Pando-derived reaction confidence, and meta-module inference before returning an artifact.
- Added a hard all-strata/all-sample barrier. Global recalibration and GEM construction are blocked when any retained stratum fails or any biological sample is absent.
- Recomputed expression capacity after the barrier by combining Human-GEM GPR-gene logCPM across all metacells, recalculating one global gene-score and GPR reaction-capacity matrix, and applying one reaction-wise Q95 calibration.
- Derived metacell-specific reaction confidence from significant Pando regions and their TF-IDF accessibility; no separate Signac `LinkPeaks()` pass is run by the integrated workflow.
- Unioned all strict-stratum reaction envelopes into one canonical `GLOBAL_UNION` meta-module and completed one shared GEM per medium scenario.
- Required condition-invariant medium constraints for shared-GEM comparison and preserved explicit no-constraint scenarios across repeated normalization.
- Changed Layer 2 parallelization to one shared-model/medium × metacell task, with one metacell-specific penalty vector and all target directions evaluated after loading the model once.
- Explicitly releases the upstream worker pool before global processing and creates a fresh Layer 2 pool; cleanup also runs on errors.

## API cleanup and contracts

- Removed the obsolete staged runner `rc_run_regcompass_multiome_metacell()`.
- Removed the sample-specific `rc_build_meta_module_gem_cache()` implementation.
- Removed `rc_load_metacell_object_from_run()` and the versioned `R/regcompass_v13.R` file.
- Removed the obsolete patch-application GitHub Actions workflow.
- Corrected parent-GEM failure semantics so an allowed reaction direction that is structurally infeasible is reported as `parent_blocked`, not `no_allowed_direction`.
- Hardened labeled microCOMPASS row-ID validation and preserved gene names in single-gene reaction-capacity calculations.
- Enforced one active pool per cell, made hard-min GPR capacity return missing for absent required subunits, restored explicit partial-GPR confidence thresholds, filtered missing pathway cross-references before meta-module expansion, and prioritized curated reaction roles in Layer 2 classification.
- Removed incomplete GPR candidates from the hard-core anchor set while retaining reactions that are independently reached through a valid hard-core subsystem or database cross-reference.
- Passes the filtered hard-core table—not the mixed candidate table—into subsystem and cross-database meta-module expansion.

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
