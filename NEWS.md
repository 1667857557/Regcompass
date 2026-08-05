# RegCompassR 2.4.0

- Defines three mutually exclusive Layer 2 structural routes: add-only FASTCORE, original MATLAB CORDA2, and a COMPASS-style full GEM.
- Applies medium scenarios as exchange-reaction bound changes only; medium handling itself never removes reaction or metabolite columns in any of the three routes.
- Retains every reaction and requested target in full-GEM mode, computes directional maximum flux under the medium, and skips the penalty LP for medium-infeasible directions.
- Records `model_completion = "none"` for the full-GEM route and rejects FASTCORE/CORDA2 controls instead of silently accepting or mislabelling them.
- Fingerprints the exact medium bounds used by full-GEM caches so same-named media with different exchange limits cannot reuse stale models.
- Uses the latest default branches of Pando_regcompass and SuperCell_Seurat_V4 without fixed revisions.
- Routes each retained cell type independently to condition GRN or standard Pando according to its retained condition count.
- Uses the current Pando condition-fit API.
- Removes retired projection and penalty fields from Layer 1, Layer 2, documentation, and result schemas.
- Consolidates mathematical definitions in `docs/mathematical-model.md`.
- Reduces user tutorials and continuous integration to the current supported workflow.
