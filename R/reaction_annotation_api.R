# Public reaction-annotation construction API.

#' Build formal reaction annotations and evidence provenance
#'
#' Constructs reaction names, formulas, direction-specific substrates and
#' products, GPR genes, database identifiers, and condition-by-cell-type
#' evidence classes. For Layer 1 schema v6, `RNA+ATAC` is assigned only when
#' Pando changes the GPR-aggregated quantitative reaction expression that is
#' actually supplied to the COMPASS-like LP penalty. Quantitative RNA is the
#' equal-weight SuperCell mean of per-cell linear CPM. Legacy Layer 1 objects
#' without quantitative reaction matrices fall back to the historical bounded
#' reaction-capacity definition.
#'
#' @param gem A validated RegCompass GEM.
#' @param layer1 Optional Layer 1 result used for GPR and evidence provenance.
#' @param reaction_ids Optional reactions to annotate.
#' @param condition_col,celltype_col Optional Layer 1 metadata columns.
#' @param evidence_tolerance Non-negative numerical comparison tolerance.
#' @return A `regcompass_reaction_annotations` list containing `reactions`,
#'   `evidence`, and `params`.
#' @export
rc_build_reaction_annotations <- function(
    gem, layer1 = NULL, reaction_ids = NULL,
    condition_col = NULL, celltype_col = NULL,
    evidence_tolerance = 1e-8) {
  if (!is.numeric(evidence_tolerance) || length(evidence_tolerance) != 1L ||
      !is.finite(evidence_tolerance) || evidence_tolerance < 0) {
    stop(
      "`evidence_tolerance` must be one non-negative finite number.",
      call. = FALSE
    )
  }
  catalog <- .rc_ra_reaction_catalog(gem, layer1, reaction_ids)
  evidence <- .rc_ra_group_evidence(
    catalog = catalog,
    layer1 = layer1,
    condition_col = condition_col,
    celltype_col = celltype_col,
    evidence_tolerance = evidence_tolerance
  )
  answer <- list(
    reactions = catalog,
    evidence = evidence,
    params = list(
      condition_col = condition_col,
      celltype_col = celltype_col,
      evidence_tolerance = evidence_tolerance,
      evidence_definition = paste(
        "RNA+ATAC requires Pando to change the GPR-aggregated quantitative",
        "reaction expression used by the Layer 2 LP penalty relative to the",
        "matched quantitative RNA-only reaction expression"
      )
    )
  )
  class(answer) <- c("regcompass_reaction_annotations", "list")
  answer
}
