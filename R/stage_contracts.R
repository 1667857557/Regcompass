.rc_stage_gem_fingerprint <- function(gem) {
  .rc_full_gem_cache_fingerprint(gem)
}

.rc_validate_metacell_artifact_contract <- function(x, argument = "metacells") {
  params <- x$params
  pooled <- x$pooled
  design <- pooled$input_design
  contract <- pooled$cache_contract
  valid <- is.list(params) && is.list(pooled) && is.list(design) &&
    is.list(contract)
  expected_label <- if (is.list(params)) {
    trimws(as.character(params$celltype_col %||% ""))
  } else {
    ""
  }
  expected_gamma <- if (is.list(params) && is.list(params$metacell_args)) {
    suppressWarnings(as.integer(params$metacell_args$gamma))
  } else {
    NA_integer_
  }
  if (valid) {
    design_label <- trimws(as.character(design$supercell_label_col %||% ""))
    design_assignment <- as.character(design$celltype_assignment %||% "")
    design_gamma <- suppressWarnings(as.integer(design$gamma))
    contract_gamma <- suppressWarnings(as.integer(
      contract$analysis_args$gamma %||% NA_integer_
    ))
    valid <- identical(
      as.character(contract$schema_version %||% ""),
      "regcompass_condition_celltype_metacell_cache_v3"
    ) &&
      isTRUE(design$condition_celltype_stratification) &&
      !isTRUE(design$condition_only_stratification) &&
      nzchar(expected_label) && !nzchar(design_label) &&
      length(design_assignment) == 1L &&
      grepl("hard condition-by-cell-type", design_assignment, fixed = TRUE) &&
      length(expected_gamma) == 1L && !is.na(expected_gamma) &&
      identical(design_gamma, expected_gamma) &&
      identical(contract_gamma, expected_gamma) &&
      isTRUE(design$depth_balance) &&
      isTRUE(contract$analysis_args$depth_balance) &&
      identical(
        as.character(design$depth_balance_policy),
        paste(
          "SuperCell local-state simplification with shared cell-type RNA UMI,",
          "ATAC fragment, and cell-count targets"
        )
      ) &&
      identical(as.character(contract$condition_col), params$condition_col) &&
      identical(as.character(contract$celltype_col), params$celltype_col) &&
      identical(as.character(contract$rna_assay), params$rna_assay) &&
      identical(as.character(contract$atac_assay), params$atac_assay) &&
      is.null(contract$label_col)
  }
  if (!isTRUE(valid)) {
    stop(
      "`", argument, "` is a legacy or incompatible metacell artifact. ",
      "Rerun `rc_regcompass_step_metacells()` with ",
      "`metacell_args = list(overwrite = TRUE)` before using it in a ",
      "current downstream stage.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_require_stage_class <- function(x, class_name, argument, producer) {
  if (!inherits(x, class_name)) {
    stop(
      "`", argument, "` must be the output of `", producer, "()`.",
      call. = FALSE
    )
  }
  if (identical(class_name, "regcompass_metacell_step")) {
    .rc_validate_metacell_artifact_contract(x, argument = argument)
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
    stop("`", argument, "` was generated from a different GEM.", call. = FALSE)
  }
  invisible(observed)
}

.rc_require_workflow_params <- function(x, expected, argument) {
  observed <- x$workflow_params %||% x$params
  if (!is.list(observed) || !identical(observed, expected)) {
    stop("`", argument, "` uses different workflow parameters.", call. = FALSE)
  }
  invisible(TRUE)
}

.rc_layer1_unit_ids <- function(layer1) {
  expression <- layer1$reaction_expression
  meta <- layer1$unit_meta
  if (!is.numeric(expression) || is.null(dim(expression)) ||
      is.null(rownames(expression)) || is.null(colnames(expression)) ||
      anyNA(rownames(expression)) || anyNA(colnames(expression)) ||
      any(!nzchar(rownames(expression))) || any(!nzchar(colnames(expression))) ||
      anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression))) {
    stop(
      "Layer 1 reaction expression must be a numeric matrix with unique reaction and unit IDs.",
      call. = FALSE
    )
  }
  if (!is.data.frame(meta)) {
    stop("Layer 1 `unit_meta` must be a data frame.", call. = FALSE)
  }
  id_col <- if ("pool_id" %in% colnames(meta)) {
    "pool_id"
  } else if ("unit_id" %in% colnames(meta)) {
    "unit_id"
  } else {
    stop("Layer 1 `unit_meta` lacks `pool_id`/`unit_id`.", call. = FALSE)
  }
  ids <- trimws(as.character(meta[[id_col]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Layer 1 unit IDs must be unique and non-empty.", call. = FALSE)
  }
  if (!identical(colnames(expression), ids)) {
    stop(
      "Layer 1 reaction-expression columns and `unit_meta` are not identically ordered.",
      call. = FALSE
    )
  }
  ids
}

.rc_validate_layer1_stage <- function(
    layer1, workflow_params = NULL, gem = NULL,
    argument = "layer1") {
  .rc_require_stage_class(
    layer1, "regcompass_layer1_step", argument,
    "rc_regcompass_step_layer1"
  )
  ids <- .rc_layer1_unit_ids(layer1)
  required_gene_matrices <- c(
    "gene_projection_common_oof", "gene_projection_condition_full_oof",
    "gene_projection_scale",
    "gene_regulatory_reliability", "gene_regulatory_reliability_available",
    "gene_regulatory_modifier_common_oof",
    "gene_regulatory_modifier_condition_full_oof",
    "gene_support_rna", "gene_support_common_oof",
    "gene_support_condition_full_oof"
  )
  missing_gene_matrices <- required_gene_matrices[
    !vapply(required_gene_matrices, function(name) {
      !is.null(layer1[[name]]) && !is.null(dim(layer1[[name]]))
    }, logical(1))
  ]
  if (length(missing_gene_matrices)) {
    stop(
      "Layer 1 is missing audited gene matrices: ",
      paste(missing_gene_matrices, collapse = ", "), call. = FALSE
    )
  }
  reference <- layer1$gene_projection_common_oof
  if (!is.numeric(reference) || is.null(rownames(reference)) ||
      is.null(colnames(reference)) || anyDuplicated(rownames(reference)) ||
      anyDuplicated(colnames(reference)) ||
      !identical(colnames(reference), ids)) {
    stop("Layer 1 gene projections require aligned unique gene and unit IDs.",
         call. = FALSE)
  }
  for (name in setdiff(
      required_gene_matrices, "gene_projection_common_oof"
  )) {
    value <- layer1[[name]]
    if (!identical(dimnames(value), dimnames(reference))) {
      stop("Layer 1 `", name, "` is not aligned with the common OOF projection.",
           call. = FALSE)
    }
  }
  if (!is.logical(layer1$gene_regulatory_reliability_available) ||
      anyNA(layer1$gene_regulatory_reliability_available)) {
    stop("Layer 1 regulatory-reliability availability must be logical.",
         call. = FALSE)
  }
  reliability <- layer1$gene_regulatory_reliability
  available <- layer1$gene_regulatory_reliability_available
  if (any(!is.finite(reliability[available])) ||
      any(reliability[available] < 0 | reliability[available] > 1 + 1e-8)) {
    stop("Layer 1 available reliability values must lie in [0, 1].",
         call. = FALSE)
  }
  projection <- layer1$gene_projection_common_oof
  projection_scale <- layer1$gene_projection_scale
  expected_modifier <- reliability * tanh(projection / projection_scale)
  expected_modifier[
    !is.finite(projection) | !is.finite(reliability)
  ] <- NA_real_
  observed_modifier <- layer1$gene_regulatory_modifier_common_oof
  comparable <- is.finite(expected_modifier) & is.finite(observed_modifier)
  if (any(is.finite(expected_modifier) != is.finite(observed_modifier)) ||
      any(abs(
        expected_modifier[comparable] - observed_modifier[comparable]
      ) > 1e-10)) {
    stop("Layer 1 regulatory modifier is not q * tanh(G / shared scale).",
         call. = FALSE)
  }
  reaction_fields <- c(
    "reaction_expression_common_oof",
    "reaction_expression_condition_full_oof",
    "reaction_expression_rna_only"
  )
  if (!all(vapply(reaction_fields, function(name) {
      value <- layer1[[name]]
      is.numeric(value) && !is.null(dim(value)) &&
        identical(dimnames(value), dimnames(layer1$reaction_expression))
    }, logical(1)))) {
    stop("Layer 1 common, condition-full, and RNA reaction matrices differ.",
         call. = FALSE)
  }
  support_fraction_fields <- c(
    "reaction_common_support_fraction",
    "reaction_condition_full_support_fraction"
  )
  if (!all(vapply(support_fraction_fields, function(name) {
      value <- layer1[[name]]
      is.numeric(value) && !is.null(dim(value)) &&
        identical(dimnames(value), dimnames(layer1$reaction_expression)) &&
        all(!is.finite(value) | (value >= 0 & value <= 1))
    }, logical(1))) ||
      !is.data.frame(layer1$zero_pattern_diagnostics) ||
      !is.data.frame(layer1$reaction_zero_support_sensitivity) ||
      !is.data.frame(layer1$reaction_link_saturation_sensitivity)) {
    stop("Layer 1 support/zero/link diagnostics are incomplete.",
         call. = FALSE)
  }
  provenance <- layer1$projection_provenance
  if (!identical(layer1$schema_version,
                 "regcompass_condition_grn_layer1_v4") ||
      !is.list(provenance) ||
      !identical(provenance$pando_schema, "pando_condition_grn_fit_v5") ||
      !identical(
        provenance$projection_origin,
        "outer_condition_stratified_cell_oof"
      ) ||
      !isTRUE(provenance$projection_used_for_penalty) ||
      !identical(provenance$full_fit_projection_used_for_penalty, FALSE) ||
      !identical(
        provenance$supercell_membership,
        "membership_table(cell_id, metacell_id)"
      )) {
    stop("Layer 1 projection provenance is absent or incompatible.",
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
  if (!is.numeric(reference) || is.null(rownames(reference)) ||
      is.null(colnames(reference)) || anyDuplicated(rownames(reference)) ||
      anyDuplicated(colnames(reference))) {
    stop("Layer 2 penalty requires unique target and unit IDs.", call. = FALSE)
  }
  for (name in setdiff(required, "penalty")) {
    value <- layer2[[name]]
    if (!identical(dimnames(value), dimnames(reference))) {
      stop("Layer 2 `", name, "` is not aligned with `penalty`.", call. = FALSE)
    }
  }
  meta <- layer2$unit_meta
  if (!is.data.frame(meta)) {
    stop("Layer 2 `unit_meta` must be a data frame.", call. = FALSE)
  }
  id_col <- if ("pool_id" %in% colnames(meta)) {
    "pool_id"
  } else if ("unit_id" %in% colnames(meta)) {
    "unit_id"
  } else {
    stop("Layer 2 `unit_meta` lacks `pool_id`/`unit_id`.", call. = FALSE)
  }
  ids <- trimws(as.character(meta[[id_col]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      !identical(colnames(reference), ids)) {
    stop("Layer 2 matrices and `unit_meta` are not identically aligned.",
         call. = FALSE)
  }
  ids
}

.rc_validate_layer2_stage <- function(
    layer2, layer1 = NULL, workflow_params = NULL, gem = NULL,
    required_mode = NULL, argument = "layer2") {
  .rc_require_stage_class(
    layer2, "regcompass_layer2_step", argument,
    "rc_regcompass_step_layer2"
  )
  ids <- .rc_layer2_unit_ids(layer2)
  comparison_matrices <- c(
    "penalty_common", "penalty_condition_full", "penalty_rna_only",
    "penalty_depth_matched_rna",
    "penalty_common_depth_interval_rna",
    "penalty_unique_increment"
  )
  if (!all(vapply(comparison_matrices, function(name) {
      value <- layer2[[name]]
      is.numeric(value) && !is.null(dim(value)) &&
        identical(dimnames(value), dimnames(layer2$penalty))
    }, logical(1))) ||
      !identical(layer2$penalty, layer2$penalty_common) ||
      !is.list(layer2$penalty_alpha_sensitivity) ||
      !length(layer2$penalty_alpha_sensitivity) ||
      !all(vapply(layer2$penalty_alpha_sensitivity, function(value) {
        is.numeric(value) && !is.null(dim(value)) &&
          identical(dimnames(value), dimnames(layer2$penalty))
      }, logical(1))) ||
      !is.list(layer2$comparison_contract) ||
      !isTRUE(layer2$comparison_contract$exact_shared_structure) ||
      !identical(
        layer2$comparison_contract$structural_model_contract,
        layer2$structural_model_contract
      ) ||
      !is.data.frame(layer2$comparison_table) ||
      !all(c(
        "reaction_id", "direction", "medium", "cell_type", "condition",
        "metacell_id", "penalty_rna_only", "penalty_common_oof",
        "penalty_condition_full_oof", "penalty_unique_increment",
        "penalty_per_target_flux", "vmax", "projection_oof_available",
        "common_support_fraction", "condition_full_support_fraction",
        "depth_sensitivity_flag", "zero_support_sensitive",
        "link_saturation_sensitive", "alpha", "inference_class",
        "comparability_class"
      ) %in% colnames(layer2$comparison_table))) {
    stop(
      "Layer 2 lacks aligned common/full/RNA reruns on one structural model.",
      call. = FALSE
    )
  }
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer2, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer2, gem, argument)
  if (!is.null(required_mode) &&
      !identical(as.character(layer2$model_mode), required_mode)) {
    stop("`", argument, "` must use `model_mode = \"", required_mode, "\"`.",
         call. = FALSE)
  }
  if (!is.null(layer1)) {
    layer1_ids <- .rc_validate_layer1_stage(
      layer1,
      workflow_params = workflow_params,
      gem = gem,
      argument = "layer1"
    )
    if (!identical(ids, layer1_ids)) {
      stop("Layer 1 and Layer 2 contain different scoring units.", call. = FALSE)
    }
  }
  invisible(ids)
}
