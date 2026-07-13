# RegCompassR 1.2.0

## Added

- Sample-specific Pando GRN inference on retained RNA+ATAC metacells.
- Metabolic target-gene selection from the intersection of single-cell RNA features and Human-GEM GPR genes.
- Export of complete and significant TF–peak–gene coefficient tables for every sample.
- Projection of significant Pando edges into sample-specific metabolic gene networks using shared-TF and direct metabolic-TF relationships.
- GRN-defined core reaction mapping and ordered meta-module expansion through core subsystems, shared KEGG/Reactome identifiers at subsystem level, and shared master-Rhea identifiers at subsystem level.
- Optional fixed-point expansion for sensitivity analysis.
- Explicit separation of biological meta-module members from transport, exchange, and other support-only reactions added for LP feasibility.
- Human-GEM v1.2 annotation preparation retaining subsystem, KEGG, Reactome, Rhea, and master-Rhea fields.
- Integrated `rc_run_regcompass_v12()` entry point.

## Compatibility

- Pins the Seurat 4-compatible Pando fork at commit `1b5f759a36630ec34d66f995906b20496a79689c`.
- Keeps the existing Layer 1 capacity/confidence calculation and all existing microCOMPASS cache strategies available.

# RegCompassR 1.1.0

- Added strict condition × sample × cell type metacell filtering gates before and after metacell construction.
- Hardened Human-GEM archive download fallback and ZIP validation.
- Updated the formal metacell README example for the current API.
