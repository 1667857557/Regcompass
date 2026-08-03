.rc_stage_gem_fingerprint <- function(gem) {
  .rc_full_gem_cache_fingerprint(gem)
}

.rc_validate_metacell_artifact_contract <- function(x, argument = "metacells") {
  params <- x$params
  pooled <- x$pooled
  design <- pooled$input_design
  contract <- pooled$cache_contract
  valid <- is.list(params) && is.list(pooled) && is.list(design) &&
    is.list(contract) &&
    identical(
      contract$schema_version,
      "regcompass_celltype_wnn_condition_joint_cache"
    ) &&
    identical(design$native_supercell_api, "SCimplify_by_graph_group") &&
    identical(design$graph_group_argument, "cell.graph.group") &&
    identical(design$condition_argument, "cell.split.condition") &&
    identical(design$graph_method, "multimodal_WNN") &&
    identical(
      design$graph_scope,
      "one_independent_WNN_graph_per_cell_type"
    ) &&
    identical(
      design$condition_scope,
      "all_conditions_joint_within_cell_type_graph"
    ) &&
    identical(
      design$membership_split_timing,
      "after_joint_WNN_graph_clustering"
    ) &&
    identical(
      design$modality_weighting,
      "adaptive_WNN_within_cell_type"
    ) &&
    identical(design$temporary_combined_stratum, FALSE) &&
    identical(contract$native_supercell_api, design$native_supercell_api) &&
    identical(contract$graph_group_argument, design$graph_group_argument) &&
    identical(contract$condition_argument, design$condition_argument) &&
    identical(contract$graph_scope, design$graph_scope) &&
    identical(contract$condition_scope, design$condition_scope) &&
    identical(
      contract$membership_split_timing,
      design$membership_split_timing
    ) &&
    identical(contract$modality_weighting, design$modality_weighting) &&
    identical(contract$condition_col, params$condition_col) &&
    identical(contract$celltype_col, params$celltype_col) &&
    identical(contract$rna_assay, params$rna_assay) &&
    identical(contract$atac_assay, params$atac_assay)
  if (!isTRUE(valid)) {
    stop(
      "`", argument,
      "` is not a cell-type-scoped, joint-condition WNN SuperCell artifact; rerun Stage 2 with overwrite=TRUE.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_require_stage_class <- function(x, class_name, argument, producer) {
  if (!inherits(x, class_name)) {
    stop("`", argument, "` must be the output of `", producer, "()`.",
         call. = FALSE)
  }
  if (identical(class_name, "regcompass_metacell_step")) {
    .rc_validate_metacell_artifact_contract(x, argument)
  }
  invisible(TRUE)
}

.rc_require_stage_gem <- function(x, gem, argument) {
  expected <- as.character(x$gem_fingerprint %||% "")
  observed <- .rc_stage_gem_fingerprint(gem)
  if (length(expected) != 1L || !nzchar(expected)) {
    stop("`", argument, "` lacks GEM provenance.", call. = FALSE)
  }
  if (!identical(expected, observed)) {
    stop("`", argument, "` was generated from a different GEM.",
         call. = FALSE)
  }
  invisible(observed)
}

.rc_require_workflow_params <- function(x, expected, argument) {
  observed <- x$workflow_params %||% x$params
  if (!is.list(observed) || !identical(observed, expected)) {
    stop("`", argument, "` uses different workflow parameters.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_layer1_unit_ids <- function(layer1) {
  expression <- layer1$reaction_expression
  meta <- layer1$unit_meta
  if (!is.numeric(expression) || is.null(dim(expression)) ||
      is.null(rownames(expression)) || is.null(colnames(expression)) ||
      anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression))) {
    stop("Layer 1 reaction expression requires unique reaction and unit IDs.",
         call. = FALSE)
  }
  if (!is.data.frame(meta)) {
    stop("Layer 1 `unit_meta` must be a data frame.", call. = FALSE)
  }
  id_col <- if ("pool_id" %in% colnames(meta)) "pool_id" else "unit_id"
  ids <- as.character(meta[[id_col]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      !identical(colnames(expression), ids)) {
    stop("Layer 1 matrices and metadata are not identically aligned.",
         call. = FALSE)
  }
  ids
}

.rc_validate_layer1_stage <- function(
    layer1, workflow_params = NULL, gem = NULL, argument = "layer1") {
  .rc_require_stage_class(
    layer1, "regcompass_layer1_step", argument, "rc_regcompass_step_layer1"
  )
  ids <- .rc_layer1_unit_ids(layer1)
  mode <- layer1$analysis_mode
  if (!mode %in% c("condition_grn", "standard_pando")) {
    stop("Layer 1 analysis mode is invalid.", call. = FALSE)
  }
  required_gene_matrices <- c(
    "gene_projection", "gene_projection_scale",
    "gene_regulatory_reliability",
    "gene_regulatory_reliability_available",
    "gene_regulatory_modifier", "gene_support_rna",
    "gene_support_multiome"
  )
  reference <- layer1$gene_projection
  if (!is.numeric(reference) || is.null(dim(reference)) ||
      !identical(colnames(reference), ids)) {
    stop("Layer 1 regulatory projection is invalid.", call. = FALSE)
  }
  for (name in required_gene_matrices) {
    value <- layer1[[name]]
    if (is.null(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(reference))) {
      stop("Layer 1 `", name, "` is missing or misaligned.",
           call. = FALSE)
    }
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
  observed <- layer1$gene_regulatory_modifier
  finite <- is.finite(expected_modifier) & is.finite(observed)
  if (any(is.finite(expected_modifier) != is.finite(observed)) ||
      any(abs(expected_modifier[finite] - observed[finite]) > 1e-10)) {
    stop("Layer 1 modifier is not reliability*tanh(projection/scale).",
         call. = FALSE)
  }
  for (name in c("reaction_expression_rna_only")) {
    value <- layer1[[name]]
    if (!is.numeric(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(layer1$reaction_expression))) {
      stop("Layer 1 reaction matrices are misaligned.", call. = FALSE)
    }
  }
  if (!is.numeric(layer1$reaction_regulatory_support_fraction) ||
      !identical(
        dimnames(layer1$reaction_regulatory_support_fraction),
        dimnames(layer1$reaction_expression)
      )) {
    stop("Layer 1 reaction regulatory support is misaligned.",
         call. = FALSE)
  }
  retired <- c(
    "reaction_expression_condition_full_oof",
    "reaction_expression_common_oof",
    "gene_projection_condition_full_oof",
    "gene_projection_common_oof",
    "gene_projection_condition_unique_oof",
    "gene_support_condition_full_oof",
    "gene_support_common_oof",
    "gene_regulatory_modifier_condition_full_oof",
    "gene_regulatory_modifier_common_oof",
    "reaction_condition_full_support_fraction",
    "reaction_common_support_fraction"
  )
  if (any(retired %in% names(layer1))) {
    stop("Layer 1 contains retired projection routes.", call. = FALSE)
  }
  provenance <- layer1$projection_provenance
  expected <- if (identical(mode, "condition_grn")) {
    list(
      origin = "paired_cell_fixed_dictionary_glm_padj_filtered",
      projection = "padj_filtered_fixed_dictionary_condition_glm",
      nonestimable =
        "coefficient_NA_and_zero_realized_penalty_contribution",
      coefficients = TRUE
    )
  } else {
    list(
      origin = "standard_pando_full_fit",
      projection = "standard_pando_full_fit",
      nonestimable = "not_applicable_standard_pando",
      coefficients = FALSE
    )
  }
  if (!identical(layer1$schema_version, "regcompass_regulatory_layer1_v4") ||
      !is.list(provenance) ||
      !identical(provenance$analysis_mode, mode) ||
      !identical(provenance$projection_origin, expected$origin) ||
      !identical(provenance$projection_name, expected$projection) ||
      !identical(provenance$nonestimable_edge_policy, expected$nonestimable) ||
      !isTRUE(provenance$projection_used_for_penalty) ||
      !identical(
        provenance$condition_coefficients_calculated,
        expected$coefficients
      ) ||
      !identical(
        provenance$supercell_membership,
        "membership_table(cell_id, metacell_id)"
      )) {
    stop("Layer 1 regulatory projection provenance is incompatible.",
         call. = FALSE)
  }
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer1, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer1, gem, argument)
  invisible(ids)
}

.rc_layer2_unit_ids <- function(layer2) {
  required <- c("penalty", "vmax", "feasible", "evaluated", "score")
  missing <- required[!vapply(required, function(name) {
    !is.null(layer2[[name]]) && !is.null(dim(layer2[[name]]))
  }, logical(1))]
  if (length(missing)) {
    stop("Layer 2 is missing matrices: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  reference <- layer2$penalty
  for (name in setdiff(required, "penalty")) {
    if (!identical(dimnames(layer2[[name]]), dimnames(reference))) {
      stop("Layer 2 matrices are not aligned.", call. = FALSE)
    }
  }
  meta <- layer2$unit_meta
  id_col <- if ("pool_id" %in% colnames(meta)) "pool_id" else "unit_id"
  ids <- as.character(meta[[id_col]])
  if (!identical(colnames(reference), ids)) {
    stop("Layer 2 matrices and metadata are not aligned.", call. = FALSE)
  }
  ids
}

.rc_validate_layer2_stage <- function(
    layer2, layer1 = NULL, workflow_params = NULL, gem = NULL,
    argument = "layer2", required_mode = NULL) {
  .rc_require_stage_class(
    layer2, "regcompass_layer2_step", argument, "rc_regcompass_step_layer2"
  )
  ids <- .rc_layer2_unit_ids(layer2)
  valid_modes <- c("meta_module_gem", "full_gem")
  if (!is.character(layer2$model_mode) || length(layer2$model_mode) != 1L ||
      is.na(layer2$model_mode) || !layer2$model_mode %in% valid_modes) {
    stop("Layer 2 model mode is invalid.", call. = FALSE)
  }
  if (!is.null(required_mode)) {
    if (!is.character(required_mode) || length(required_mode) != 1L ||
        is.na(required_mode) || !required_mode %in% valid_modes) {
      stop(
        "`required_mode` must be one of: ",
        paste(valid_modes, collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (!identical(layer2$model_mode, required_mode)) {
      stop(
        "Layer 2 model mode is `", layer2$model_mode,
        "`; `", required_mode, "` is required.", call. = FALSE
      )
    }
  }
  reference <- layer2$penalty
  for (name in c("penalty_rna_only", "score_rna_only")) {
    value <- layer2[[name]]
    if (!is.numeric(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(reference))) {
      stop("Layer 2 route `", name, "` is missing or misaligned.",
           call. = FALSE)
    }
  }
  retired <- c(
    "penalty_condition_full_oof", "penalty_common_oof",
    "penalty_condition_unique_increment",
    "score_condition_full_oof_display_only",
    "score_common_oof_display_only",
    "score_rna_only_display_only"
  )
  if (any(retired %in% names(layer2))) {
    stop("Layer 2 contains retired penalty routes.", call. = FALSE)
  }
  contract <- layer2$comparison_contract
  required_contract <- c(
    "primary", "rna_control", "nonestimable_edge_policy",
    "exact_shared_structure", "structural_model_contract",
    "effect_size_basis", "ecdf_effect_size_eligible"
  )
  if (!identical(layer2$schema_version, "regcompass_regulatory_layer2_v3") ||
      !is.list(contract) ||
      !all(required_contract %in% names(contract)) ||
      !identical(contract$primary, "penalty") ||
      !identical(contract$rna_control, "penalty_rna_only") ||
      !isTRUE(contract$exact_shared_structure) ||
      !identical(
        contract$nonestimable_edge_policy,
        "coefficient_NA_and_zero_realized_penalty_contribution"
      )) {
    stop("Layer 2 comparison contract is incomplete.", call. = FALSE)
  }
  if (!is.null(layer1)) {
    .rc_validate_layer1_stage(layer1, argument = "layer1")
    if (!identical(ids, .rc_layer1_unit_ids(layer1))) {
      stop("Layer 1 and Layer 2 unit IDs differ.", call. = FALSE)
    }
  }
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer2, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer2, gem, argument)
  invisible(ids)
}
