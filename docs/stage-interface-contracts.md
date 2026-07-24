# Stage input-output contracts

RegCompassR 1.8.3 connects stages only when their classes, workflow settings, GEM provenance, and scoring-unit order agree.

| Function | Required input | Output class | Downstream contract |
|---|---|---|---|
| `rc_regcompass_step_grn()` | paired RNA+ATAC object, GEM, motifs, genome | `regcompass_grn_step` | stores condition/cell-type/assay settings and GEM fingerprint |
| `rc_regcompass_step_metacells()` | original paired object and the same metadata-column names | `regcompass_metacell_step` | stores metacell object, membership, audited labels, and workflow parameters |
| `rc_regcompass_step_meta_modules()` | matching GRN and metacell stages plus the Stage 1 GEM | `regcompass_meta_module_step` | verifies group coverage and GEM fingerprint; stores global core and union memberships |
| `rc_regcompass_step_layer1()` | metacell and meta-module stages from the same workflow and GEM | `regcompass_layer1_step` | reaction-expression columns must be identical to ordered `unit_meta$pool_id` |
| `rc_regcompass_step_layer2()` | Layer 1, matching global modules, GEM, and shared medium | `regcompass_layer2_step` | matrices share identical target/unit dimnames; stores core set, workflow parameters, GEM fingerprint, and persistent model files |
| `rc_regcompass_step_target_union()` | matching Stage 3-5 objects from a union-GEM run | `regcompass_target_union_step` | anchors must be previous core targets; second-pass targets must be non-core reactions directly sharing KEGG, Reactome, or master-Rhea IDs with an anchor and must exist in the original cached union GEM |
| `rc_regcompass_step_results()` | matching Stage 1-5 objects and GEM | final result list | rejects different GEMs, workflow settings, classes, or unit order before ranking and annotation |

## Files that must persist

- Each stage wrapper RDS used for restart.
- Stage 5 files listed in `model_cache_summary$file`.
- The same GEM object or a GEM with an identical fingerprint.

Compact inspection artifacts do not replace their stage wrapper because they do not carry the full input contract.

## Fail-fast conditions

The workflow stops when:

- a required stage class is absent;
- condition, cell-type, or assay settings differ;
- the GEM fingerprint differs;
- Layer 1 and Layer 2 scoring units are missing, duplicated, or reordered;
- the Stage 5 core set differs from Stage 3;
- selected cores have no direct KEGG, Reactome, or master-Rhea-linked non-core reactions;
- a directly mapped target is absent from the original union;
- a cached union-model file is missing or has changed provenance.

Same-subsystem, fixed-point, and transitive target-union expansion are not supported. These checks prevent numerically valid but biologically unrelated stage objects from being combined.
