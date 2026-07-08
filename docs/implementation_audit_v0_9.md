# Documentation audit

Checked `README.md` against the current R API.

## Public functions used in the tutorial

- `rc_validate_multiome_input()` — current arguments match the example.
- `rc_prepare_human2_gem()` — example uses pinned `version`; `save_rds` is optional.
- `rc_annotate_reaction_roles()` — example uses optional `reaction_role_table`.
- `rc_make_medium_scenarios()` — example uses supported `blood_like` scenario.
- `rc_run_regcompass_multiome_metacell()` — example uses current metacell workflow arguments.
- `rc_select_target_reactions()` — example uses current selection arguments.
- `rc_run_microcompass()` — example uses current strict microCOMPASS arguments.
- `rc_test_microcompass_differential()` — example uses current differential-test arguments.
- `rc_export_microcompass()` — example uses current export signature.

## Cleanup result

The tutorial now lists only the current strict multiome workflow and exported functions used by that workflow. Old aliases, abandoned analysis branches, and internal helper APIs are omitted.
