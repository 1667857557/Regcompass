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
      "regcompass_shared_walktrap_condition_cut_cache_v1"
    ) &&
    identical(contract$native_supercell_api, "SCimplify_by_graph_group") &&
    identical(contract$graph_group_argument, "cell.graph.group") &&
    identical(contract$condition_argument, "cell.split.condition") &&
    identical(contract$condition_partition, "hierarchy_constrained") &&
    identical(
      contract$partition_schema_version,
      "shared_walktrap_condition_cut_v1"
    ) &&
    identical(
      contract$graph_scope,
      "one_independent_WNN_graph_per_cell_type"
    ) &&
    identical(
      contract$condition_scope,
      "shared_WNN_and_Walktrap_with_condition_specific_hierarchy_cut"
    ) &&
    identical(
      contract$membership_split_timing,
      "condition_specific_cut_of_shared_walktrap_hierarchy"
    ) &&
    identical(
      contract$graph_method,
      "SuperCell_multimodal_WNN_then_walktrap"
    ) &&
    identical(
      contract$aggregation_method,
      "SCimplify_for_Seurat_membership_mode"
    ) &&
    identical(
      contract$modality_weighting,
      "adaptive_WNN_within_cell_type"
    ) &&
    identical(design$native_supercell_api, "SCimplify_by_graph_group") &&
    identical(design$graph_group_argument, "cell.graph.group") &&
    identical(design$condition_argument, "cell.split.condition") &&
    identical(design$condition_partition, "hierarchy_constrained") &&
    identical(
      design$partition_schema_version,
      "shared_walktrap_condition_cut_v1"
    ) &&
    identical(design$graph_method, "multimodal_WNN") &&
    identical(
      design$clustering_method,
      "one_shared_walktrap_hierarchy_per_cell_type"
    ) &&
    identical(
      design$final_partition_method,
      "condition_specific_finest_feasible_cut_of_shared_hierarchy"
    ) &&
    identical(
      design$aggregation_method,
      "SCimplify_for_Seurat_with_membership"
    ) &&
    identical(
      design$graph_scope,
      "one_independent_WNN_graph_per_cell_type"
    ) &&
    identical(
      design$condition_scope,
      "all_conditions_joint_for_WNN_and_Walktrap_then_condition_specific_hierarchy_cut"
    ) &&
    identical(
      design$membership_split_timing,
      "during_final_shared_hierarchy_cut_selection"
    ) &&
    identical(
      design$modality_weighting,
      "adaptive_WNN_within_cell_type"
    ) &&
    identical(design$temporary_combined_stratum, FALSE) &&
    identical(pooled$partition_policy, "hierarchy_constrained") &&
    identical(
      pooled$partition_schema_version,
      design$partition_schema_version
    ) &&
    identical(
      contract$partition_schema_version,
      pooled$partition_schema_version
    ) &&
    identical(contract$condition_col, params$condition_col) &&
    identical(contract$celltype_col, params$celltype_col) &&
    identical(contract$rna_assay, params$rna_assay) &&
    identical(contract$atac_assay, params$atac_assay)
  if (!isTRUE(valid)) {
    stop(
      "`", argument,
      "` is not a shared-Walktrap hierarchy-constrained SuperCell artifact; rerun Stage 2 with the current RegCompass/SuperCell contract.",
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
  valid_modes <- c("condition_grn", "standard_pando", "mixed_pando")
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% valid_modes) {
    stop("Layer 1 analysis mode is invalid.", call. = FALSE)
  }

  required_gene_matrices <- c(
    "gene_projection", "gene_projection_scale",
    "gene_regulatory_reliability",
    "gene_regulatory_reliability_available",
    "gene_regulatory_modifier", "gene_support_rna",
    "gene_support_multiome", "gene_expression_quantitative_rna",
    "gene_expression_quantitative_multiome",
    "rna_metacell_mean_single_cell_cpm"
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
  if (!identical(
        layer1$rna_metacell_mean_single_cell_cpm,
        layer1$gene_expression_quantitative_rna
      )) {
    stop(
      "Layer 1 quantitative RNA must equal the saved SuperCell mean single-cell CPM matrix.",
      call. = FALSE
    )
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
    stop("Layer 1 regulatory modifier is inconsistent with its inputs.",
         call. = FALSE)
  }

  quantitative_rna <- layer1$gene_expression_quantitative_rna
  if (any(quantitative_rna[is.finite(quantitative_rna)] < 0)) {
    stop("Layer 1 quantitative RNA expression must be non-negative.",
         call. = FALSE)
  }
  expected_quantitative <- .rc_integrate_regulatory_expression(
    quantitative_rna, observed
  )
  observed_quantitative <- layer1$gene_expression_quantitative_multiome
  finite <- is.finite(expected_quantitative) & is.finite(observed_quantitative)
  if (any(is.finite(expected_quantitative) !=
          is.finite(observed_quantitative)) ||
      any(abs(expected_quantitative[finite] - observed_quantitative[finite]) >
          1e-10 * pmax(1, abs(expected_quantitative[finite])))) {
    stop(
      "Layer 1 quantitative multiome expression is inconsistent with X*2^R.",
      call. = FALSE
    )
  }

  required_reaction_matrices <- c(
    "reaction_expression_rna_only",
    "reaction_expression_quantitative",
    "reaction_expression_quantitative_rna_only",
    "reaction_structural_support",
    "reaction_structural_support_rna_only"
  )
  for (name in required_reaction_matrices) {
    value <- layer1[[name]]
    if (!is.numeric(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(layer1$reaction_expression))) {
      stop("Layer 1 `", name, "` is missing or misaligned.",
           call. = FALSE)
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
  if (!identical(
        layer1$reaction_expression,
        layer1$reaction_structural_support
      ) ||
      !identical(
        layer1$reaction_expression_rna_only,
        layer1$reaction_structural_support_rna_only
      )) {
    stop(
      "Layer 1 compatibility reaction_expression fields must remain bounded structural support.",
      call. = FALSE
    )
  }

  route_attr <- "regcompass_quantitative_penalty_route"
  if (!identical(attr(layer1$reaction_expression, route_attr, exact = TRUE),
                 "multiome") ||
      !identical(attr(layer1$reaction_expression_rna_only,
                      route_attr, exact = TRUE), "rna_only")) {
    stop(
      "Layer 1 structural compatibility matrices lack deterministic quantitative penalty route markers.",
      call. = FALSE
    )
  }
  if (!is.logical(layer1$reaction_expression_available) ||
      !identical(
        dimnames(layer1$reaction_expression_available),
        dimnames(layer1$reaction_expression)
      ) ||
      !identical(
        unname(layer1$reaction_expression_available),
        unname(is.finite(layer1$reaction_expression))
      ) ||
      !is.logical(layer1$reaction_expression_quantitative_available) ||
      !identical(
        dimnames(layer1$reaction_expression_quantitative_available),
        dimnames(layer1$reaction_expression_quantitative)
      ) ||
      !identical(
        unname(layer1$reaction_expression_quantitative_available),
        unname(is.finite(layer1$reaction_expression_quantitative))
      )) {
    stop("Layer 1 reaction-expression availability masks are inconsistent.",
         call. = FALSE)
  }

  quantitative_contract <- layer1$quantitative_penalty_contract
  structural_contract <- layer1$structural_support_contract
  if (!is.list(quantitative_contract) ||
      !isTRUE(quantitative_contract$bounded_support_excluded_from_lp_penalty) ||
      !identical(
        quantitative_contract$baseline_gene_expression,
        "equal_mean_single_cell_linear_cpm"
      ) ||
      !identical(
        quantitative_contract$single_cell_normalization,
        "raw_counts_per_cell_divided_by_complete_cell_RNA_library_times_1e6"
      ) ||
      !identical(
        quantitative_contract$metacell_aggregation,
        "equal_mean_by_exact_SuperCell_membership"
      ) ||
      !identical(
        quantitative_contract$library_size_weighted_metacell_average,
        FALSE
      ) ||
      !identical(quantitative_contract$regulatory_multiplier, "2^R") ||
      !identical(quantitative_contract$structural_route_marker, route_attr) ||
      !identical(quantitative_contract$primary_route, "multiome") ||
      !identical(quantitative_contract$rna_control_route, "rna_only") ||
      !is.list(structural_contract) ||
      !identical(
        structural_contract$intended_use,
        "CORDA2_and_structural_confidence"
      ) ||
      !identical(structural_contract$quantitative_lp_penalty, FALSE) ||
      !isTRUE(structural_contract$latent_cpm_structural_only)) {
    stop("Layer 1 quantitative/structural support contracts are incomplete.",
         call. = FALSE)
  }

  provenance <- layer1$projection_provenance
  required_provenance <- c(
    "analysis_mode", "pando_schema", "projection_origin",
    "projection_used_for_penalty", "projection_name",
    "condition_coefficients_calculated", "supercell_membership",
    "quantitative_rna_source", "quantitative_rna_aggregation",
    "unavailable_target_policy", "nonestimable_edge_policy",
    "cell_type_analysis_mode"
  )
  routing <- provenance$cell_type_analysis_mode
  if (!identical(layer1$schema_version, "regcompass_regulatory_layer1_v6") ||
      !is.list(provenance) ||
      !all(required_provenance %in% names(provenance)) ||
      !identical(provenance$analysis_mode, mode) ||
      !isTRUE(provenance$projection_used_for_penalty) ||
      !identical(
        provenance$supercell_membership,
        "membership_table(cell_id, metacell_id)"
      ) ||
      !identical(
        provenance$quantitative_rna_source,
        "Pando_retained_cell_level_raw_RNA_counts"
      ) ||
      !identical(
        provenance$quantitative_rna_aggregation,
        "equal_mean_after_per_cell_linear_CPM_normalization"
      ) ||
      !is.data.frame(routing) ||
      !all(c("cell_type", "analysis_mode") %in% colnames(routing)) ||
      !nrow(routing) ||
      anyNA(routing$cell_type) || anyNA(routing$analysis_mode) ||
      any(!routing$analysis_mode %in% c("condition_grn", "standard_pando")) ||
      any(!nzchar(as.character(provenance$pando_schema))) ||
      any(!nzchar(as.character(provenance$projection_origin))) ||
      any(!nzchar(as.character(provenance$projection_name))) ||
      any(!nzchar(as.character(provenance$nonestimable_edge_policy)))) {
    stop("Layer 1 regulatory projection provenance is incompatible.",
         call. = FALSE)
  }
  observed_modes <- sort(unique(as.character(routing$analysis_mode)))
  expected_modes <- switch(
    mode,
    condition_grn = "condition_grn",
    standard_pando = "standard_pando",
    mixed_pando = sort(c("condition_grn", "standard_pando"))
  )
  if (!identical(observed_modes, expected_modes)) {
    stop("Layer 1 cell-type routing does not match its analysis mode.",
         call. = FALSE)
  }
  expected_coefficients <- "condition_grn" %in% observed_modes
  if (!identical(
        provenance$condition_coefficients_calculated,
        expected_coefficients
      )) {
    stop("Layer 1 condition-coefficient provenance is inconsistent.",
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

.rc_validate_layer2_stage_core <- function(
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
      stop("`required_mode` must be one of: ",
           paste(valid_modes, collapse = ", "), ".", call. = FALSE)
    }
    if (!identical(layer2$model_mode, required_mode)) {
      stop("Layer 2 model mode is `", layer2$model_mode,
           "`; `", required_mode, "` is required.", call. = FALSE)
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
      !identical(layer2$evidence_policy, "quantitative_penalty_only") ||
      !identical(
        contract$nonestimable_edge_policy,
        "coefficient_NA_and_zero_realized_penalty_contribution"
      )) {
    stop(
      "Layer 2 comparison/quantitative penalty contract is incomplete; rerun Stage 5 with the current RegCompass version.",
      call. = FALSE
    )
  }
  if (!is.null(layer1)) {
    .rc_validate_layer1_stage(layer1, argument = "layer1")
    if (!identical(ids, .rc_layer1_unit_ids(layer1))) {
      stop("Layer 1 and Layer 2 unit IDs differ.", call. = FALSE)
    }

    components <- layer2$penalty_components$reaction_expression
    expected <- layer1$reaction_expression_quantitative
    if (!is.numeric(components) || is.null(dim(components)) ||
        !identical(colnames(components), ids)) {
      stop("Layer 2 quantitative reaction-expression provenance is missing.",
           call. = FALSE)
    }
    common <- intersect(rownames(expected), rownames(components))
    if (!length(common)) {
      stop("Layer 1 and Layer 2 share no quantitative reaction-expression rows.",
           call. = FALSE)
    }
    observed_input <- components[common, ids, drop = FALSE]
    expected_input <- expected[common, ids, drop = FALSE]
    finite <- is.finite(observed_input) & is.finite(expected_input)
    if (any(is.finite(observed_input) != is.finite(expected_input)) ||
        any(abs(observed_input[finite] - expected_input[finite]) >
            1e-10 * pmax(1, abs(expected_input[finite])))) {
      stop(
        "Layer 2 primary penalty was not constructed from Layer 1 quantitative multiome reaction expression; rerun Stage 5.",
        call. = FALSE
      )
    }
  }
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer2, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer2, gem, argument)
  invisible(ids)
}

# Progress-aware entry point; the algorithm remains in the core above.
.rc_validate_layer2_stage <- function(...) {
  answer <- do.call(
    .rc_validate_layer2_stage_core,
    list(...)
  )
  .rc_layer2_overall_event(
    "layer2_validation_complete", 10L,
    "Layer 2 schemas, ordering and shared-model contracts validated"
  )
  answer
}
