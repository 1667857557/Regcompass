# RegCompassR documentation audit

The user-facing tutorial is `README.md`.

## Current public workflow

1. Validate multiome input with `rc_validate_multiome_input()`.
2. Build or import metacells with `rc_make_metacells()` or `rc_import_metacells()`.
3. Read a provenance-tracked GEM with `rc_read_gem()` and annotate roles with `rc_annotate_reaction_roles()`.
4. Create explicit medium scenarios with `rc_make_medium_scenarios()`.
5. Run Layer 1 with `rc_run_layer1_multiome()`.
6. Select targets with `rc_select_target_reactions()`.
7. Run cached structural micro-GEM strict LP with `rc_run_microcompass()`.
8. Test and export with `rc_test_microcompass_differential()` and `rc_export_microcompass()`.

## Removed from the tutorial

- Long historical implementation plan text.
- Pool/pseudobulk explanations that are not the main workflow.
- Claims that metacells are biological replicates.
- Relaxed LP/FVA tutorial steps.
- Full-GEM flux, uptake/secretion, enzyme activity, or ATAC-causality wording.
