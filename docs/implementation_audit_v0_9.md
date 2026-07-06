# RegCompassR documentation audit

This file is intentionally short. The user-facing tutorial is `README.md`.

## Current public workflow

1. `rc_validate_multiome_input()` checks annotated Seurat/Signac raw-count input.
2. `rc_make_metacells()` or `rc_import_metacells()` supplies metacell counts and metadata.
3. `rc_run_layer1_multiome()` computes Layer 1 reaction evidence.
4. `rc_read_gem()`, `rc_validate_gem()`, `rc_annotate_reaction_roles()`, and `rc_apply_medium_constraints()` prepare the GEM.
5. `rc_select_target_reactions()` chooses target reactions only.
6. `rc_build_target_microgem()` builds target-local networks.
7. `rc_run_microcompass()` runs strict LP and optional relaxed LP/FVA.
8. `rc_test_microcompass_differential()` and `rc_export_microcompass()` summarize and export results.

## Removed from the tutorial

- Long historical pool/pseudobulk explanations.
- Claims that metacells are biological replicates.
- Full-GEM flux, uptake/secretion, enzyme activity, or causality wording.
- Detailed implementation claims not exposed by current function signatures.
