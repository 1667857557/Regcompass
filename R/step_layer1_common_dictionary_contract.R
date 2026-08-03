# Layer 1 validation for the common-dictionary condition projection.

.rc_validate_layer1_stage_original <- .rc_validate_layer1_stage

.rc_validate_layer1_stage <- function(
    layer1, workflow_params = NULL, gem = NULL, argument = "layer1") {
  mode <- layer1$analysis_mode %||% NA_character_
  if (!identical(mode, "condition_grn")) {
    return(.rc_validate_layer1_stage_original(
      layer1 = layer1,
      workflow_params = workflow_params,
      gem = gem,
      argument = argument
    ))
  }
  .rc_require_stage_class(
    layer1, "regcompass_layer1_step", argument, "rc_regcompass_step_layer1"
  )
  ids <- .rc_layer1_unit_ids(layer1)
  required_gene_matrices <- c(
    "gene_projection_condition_full_oof",
    "gene_projection_common_oof",
    "gene_projection_condition_unique_oof",
    "gene_projection_scale",
    "gene_regulatory_reliability",
    "gene_regulatory_reliability_available",
    "gene_regulatory_modifier_condition_full_oof",
    "gene_regulatory_modifier_common_oof",
    "gene_support_rna",
    "gene_support_condition_full_oof",
    "gene_support_common_oof"
  )
  reference <- layer1$gene_projection_condition_full_oof
  if (!is.numeric(reference) || is.null(dim(reference)) ||
      !identical(colnames(reference), ids)) {
    stop("Layer 1 primary common-dictionary projection is invalid.",
         call. = FALSE)
  }
  for (name in required_gene_matrices) {
    value <- layer1[[name]]
    if (is.null(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(reference))) {
      stop("Layer 1 `", name, "` is missing or misaligned.",
           call. = FALSE)
    }
  }
  if (!isTRUE(all.equal(
        layer1$gene_projection_common_oof,
        reference,
        tolerance = 1e-12,
        check.attributes = TRUE
      )) ||
      !isTRUE(all.equal(
        layer1$gene_projection_condition_unique_oof,
        reference * 0,
        tolerance = 1e-12,
        check.attributes = TRUE
      ))) {
    stop("Layer 1 compatibility projection aliases are inconsistent.",
         call. = FALSE)
  }
  if (!is.logical(layer1$gene_regulatory_reliability_available) ||
      anyNA(layer1$gene_regulatory_reliability_available)) {
    stop("Layer 1 reliability availability must be logical.",
         call. = FALSE)
  }
  expected_modifier <- layer1$gene_regulatory_reliability *
    tanh(reference / layer1$gene_projection_scale)
  expected_modifier[
    !is.finite(reference) |
      !is.finite(layer1$gene_regulatory_reliability)
  ] <- NA_real_
  observed <- layer1$gene_regulatory_modifier_condition_full_oof
  finite <- is.finite(expected_modifier) & is.finite(observed)
  if (any(is.finite(expected_modifier) != is.finite(observed)) ||
      any(abs(expected_modifier[finite] - observed[finite]) > 1e-10)) {
    stop("Layer 1 modifier is not reliability*tanh(projection/scale).",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(
        layer1$gene_regulatory_modifier_common_oof,
        observed,
        tolerance = 1e-12,
        check.attributes = TRUE
      ))) {
    stop("Layer 1 common modifier compatibility alias is inconsistent.",
         call. = FALSE)
  }
  for (name in c(
    "reaction_expression_condition_full_oof",
    "reaction_expression_common_oof",
    "reaction_expression_rna_only"
  )) {
    value <- layer1[[name]]
    if (!is.numeric(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(layer1$reaction_expression))) {
      stop("Layer 1 reaction matrices are misaligned.", call. = FALSE)
    }
  }
  if (!identical(
        layer1$reaction_expression,
        layer1$reaction_expression_condition_full_oof
      ) ||
      !identical(
        layer1$reaction_expression_common_oof,
        layer1$reaction_expression_condition_full_oof
      )) {
    stop("Layer 1 reaction compatibility aliases are inconsistent.",
         call. = FALSE)
  }
  provenance <- layer1$projection_provenance
  if (!identical(layer1$schema_version, "regcompass_regulatory_layer1_v3") ||
      !is.list(provenance) ||
      !identical(provenance$analysis_mode, "condition_grn") ||
      !identical(
        provenance$projection_origin,
        "paired_cell_full_fit_fixed_dictionary_glm_padj_filtered"
      ) ||
      !identical(
        provenance$primary_projection,
        "padj_filtered_fixed_dictionary_condition_glm"
      ) ||
      !identical(
        provenance$nonestimable_edge_policy,
        "coefficient_NA_and_zero_realized_penalty_contribution"
      ) ||
      !isTRUE(provenance$projection_used_for_penalty) ||
      !isTRUE(provenance$condition_coefficients_calculated) ||
      !identical(
        provenance$supercell_membership,
        "membership_table(cell_id, metacell_id)"
      )) {
    stop("Layer 1 common-dictionary projection provenance is incompatible.",
         call. = FALSE)
  }
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer1, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer1, gem, argument)
  invisible(ids)
}
