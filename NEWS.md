# RegCompassR 2.4.0

- Defines three mutually exclusive Layer 2 structural builders: add-only FASTCORE, original MATLAB CORDA2, and a medium-constrained full GEM pruned only for flux consistency.
- Records `model_completion = "none"` for the full-GEM route and rejects FASTCORE/CORDA2 controls instead of silently accepting or mislabelling them.
- Uses solver- and threshold-aware cache fingerprints for medium-pruned full GEMs and exports the number of removed flux-inconsistent reactions.
- Uses the latest default branches of Pando_regcompass and SuperCell_Seurat_V4 without fixed revisions.
- Routes each retained cell type independently to condition GRN or standard Pando according to its retained condition count.
- Uses the current Pando condition-fit API.
- Removes retired projection and penalty fields from Layer 1, Layer 2, documentation, and result schemas.
- Consolidates mathematical definitions in `docs/mathematical-model.md`.
- Reduces user tutorials and continuous integration to the current supported workflow.