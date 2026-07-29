.rc_cell_first_projection_layer1 <- function(
    grn_result, metacell_object, membership, metacell_meta, gem,
    sample_col, condition_col, celltype_col, rna_assay,
    projection_component, comparison_support, regulatory_alpha,
    gpr_and_method, gene_half_saturation, parallel, BPPARAM) {
  if (!identical(
    grn_result$schema_version, "regcompass_condition_grn_fit_v4"
  ) || !inherits(grn_result$pando_grn_data, "GRNData")) {
    stop("Stage 1 lacks the exact Pando ConditionGRNFit v4 GRNData artifact.",
         call. = FALSE)
  }
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership)) ||
      anyDuplicated(as.character(membership$cell_id))) {
    stop("SuperCell membership must map each cell once using cell_id/metacell_id.",
         call. = FALSE)
  }
  units <- as.character(colnames(metacell_object))
  if (!setequal(unique(as.character(membership$metacell_id)), units)) {
    stop("SuperCell membership metacell IDs do not match the metacell object.",
         call. = FALSE)
  }
  if (!is.data.frame(metacell_meta) ||
      !all(c(condition_col, celltype_col) %in% colnames(metacell_meta))) {
    stop("SuperCell metadata lack condition or broad-cell-type strata.",
         call. = FALSE)
  }
  id_col <- if ("metacell_id" %in% colnames(metacell_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(metacell_meta)) {
    "pool_id"
  } else {
    stop("SuperCell metadata lack metacell_id/pool_id.", call. = FALSE)
  }
  unit_meta <- metacell_meta
  unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  if (anyNA(unit_meta$pool_id) || any(!nzchar(unit_meta$pool_id)) ||
      anyDuplicated(unit_meta$pool_id) ||
      !setequal(unit_meta$pool_id, units)) {
    stop("SuperCell metadata IDs do not match the assay columns.",
         call. = FALSE)
  }
  unit_meta <- unit_meta[
    match(units, as.character(unit_meta$pool_id)), , drop = FALSE
  ]
  unit_meta$unit_id <- unit_meta$pool_id
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  counts <- .rc_get_assay_counts(metacell_object, rna_assay)
  if (!identical(colnames(counts), units)) {
    stop("SuperCell RNA assay columns are not aligned to metacell IDs.",
         call. = FALSE)
  }
  full_library_size <- Matrix::colSums(counts)
  rna_counts <- counts[
    tolower(rownames(counts)) %in% gpr_genes, , drop = FALSE
  ]
  rna_logcpm <- .rc_metacell_logcpm(
    rna_counts, library_size = full_library_size[colnames(rna_counts)]
  )
  rownames(rna_logcpm) <- tolower(rownames(rna_logcpm))
  if (anyDuplicated(rownames(rna_logcpm))) {
    stop("Duplicated GPR genes after case normalization.", call. = FALSE)
  }
  genes <- rownames(rna_logcpm)
  empty_projection <- matrix(
    NA_real_, nrow = length(genes), ncol = length(units),
    dimnames = list(genes, units)
  )
  projection_absolute <- empty_projection
  projection_shared <- empty_projection
  projection_deviation <- empty_projection
  reliability <- empty_projection
  reliability_available <- matrix(
    FALSE, nrow = length(genes), ncol = length(units),
    dimnames = list(genes, units)
  )
  coverage <- list()
  fits <- grn_result$condition_grn_fits
  if (!is.list(fits) || !length(fits)) {
    stop("Stage 1 contains no Pando ConditionGRNFit v4 contracts.",
         call. = FALSE)
  }
  assign_projection <- function(destination, projected) {
    score <- t(as.matrix(projected$gene_score))
    normalized <- tolower(rownames(score))
    if (anyDuplicated(normalized)) {
      stop("Pando projection targets duplicate after case normalization.",
           call. = FALSE)
    }
    rownames(score) <- normalized
    target <- intersect(rownames(score), rownames(destination))
    group <- intersect(colnames(score), colnames(destination))
    if (!setequal(colnames(score), group)) {
      stop("Pando group projections contain unknown SuperCell IDs.",
           call. = FALSE)
    }
    destination[target, group] <- score[target, group, drop = FALSE]
    destination
  }
  for (fit in fits) {
    if (!identical(fit$schema_version, "pando_condition_grn_fit_v4")) {
      stop("Layer 1 requires only ConditionGRNFit v4 contracts.", call. = FALSE)
    }
    fit_comparison_support <- if (identical(comparison_support, "auto")) {
      if (length(fit$condition_levels) == 2L) {
        "pairwise_common"
      } else {
        "global_common"
      }
    } else {
      comparison_support
    }
    if (identical(fit_comparison_support, "pairwise_common") &&
        length(fit$condition_levels) != 2L) {
      stop(
        "pairwise_common requires exactly two conditions per cell type; ",
        "use global_common for a multi-condition omnibus analysis.",
        call. = FALSE
      )
    }
    project_one <- function(component) {
      Pando::project_condition_grn_groups(
        object = grn_result$pando_grn_data,
        fit = fit,
        membership = membership,
        group_col = "metacell_id",
        component = component,
        scale = "std",
        support_policy = fit_comparison_support
      )
    }
    absolute <- project_one("condition")
    shared <- project_one("shared")
    deviation <- project_one("deviation")
    projection_absolute <- assign_projection(
      projection_absolute, absolute
    )
    projection_shared <- assign_projection(
      projection_shared, shared
    )
    projection_deviation <- assign_projection(
      projection_deviation, deviation
    )
    fit_units <- rownames(absolute$gene_score)
    expected_units <- as.character(unit_meta$pool_id[
      as.character(unit_meta[[celltype_col]]) == fit$cell_type &
        as.character(unit_meta[[condition_col]]) %in%
          fit$condition_levels
    ])
    if (!setequal(fit_units, expected_units)) {
      stop(
        "Pando paired cells do not cover every SuperCell in cell type `",
        fit$cell_type, "`.", call. = FALSE
      )
    }
    pooled_oof <- fit$target_rsq_oof_pooled
    if (!is.numeric(pooled_oof) || anyDuplicated(tolower(names(pooled_oof))) ||
        is.null(names(pooled_oof)) || any(!nzchar(names(pooled_oof)))) {
      stop("ConditionGRNFit v4 pooled OOF metrics are invalid.",
           call. = FALSE)
    }
    names(pooled_oof) <- tolower(names(pooled_oof))
    blocked_available <- fit$sample_blocked_oof_available
    if (!is.logical(blocked_available) || anyNA(blocked_available) ||
        is.null(names(blocked_available)) ||
        !setequal(names(blocked_available), names(
          fit$target_rsq_oof_pooled
        ))) {
      stop(
        "ConditionGRNFit v4 sample-blocked OOF availability is invalid.",
        call. = FALSE
      )
    }
    names(blocked_available) <- tolower(names(blocked_available))
    response_independent_candidates <- identical(
      fit$candidate_screen, "motif_domain"
    )
    status <- absolute$source_projection$target_condition_status
    status_targets <- tolower(as.character(status$target))
    if (anyNA(status_targets) || any(!nzchar(status_targets)) ||
        any(!status_targets %in% names(pooled_oof)) ||
        any(!status_targets %in% names(blocked_available))) {
      stop(
        "Pando target status is not aligned to pooled OOF reliability fields.",
        call. = FALSE
      )
    }
    fit_targets <- intersect(names(pooled_oof), genes)
    q <- sqrt(pmax(0, as.numeric(pooled_oof[fit_targets])))
    available <- response_independent_candidates &
      as.logical(blocked_available[fit_targets]) &
      is.finite(as.numeric(pooled_oof[fit_targets]))
    q[!available | !is.finite(q)] <- 0
    reliability[fit_targets, fit_units] <- q
    reliability_available[fit_targets, fit_units] <- available
    status$cell_type <- fit$cell_type
    status$comparison_support <- fit_comparison_support
    status$response_independent_candidate_graph <-
      response_independent_candidates
    status$sample_blocked_oof_available <- as.logical(
      blocked_available[tolower(status$target)]
    )
    status$gene_regulatory_reliability_available <-
      status$response_independent_candidate_graph &
      status$sample_blocked_oof_available &
      is.finite(as.numeric(pooled_oof[tolower(status$target)]))
    status$regulatory_reliability <- sqrt(pmax(
      0, as.numeric(pooled_oof[tolower(status$target)])
    ))
    status$regulatory_reliability[
      !status$gene_regulatory_reliability_available |
        !is.finite(status$regulatory_reliability)
    ] <- 0
    coverage[[length(coverage) + 1L]] <- status
  }
  projection <- switch(
    projection_component,
    condition = projection_absolute,
    shared = projection_shared,
    deviation = projection_deviation
  )
  modifier <- reliability * tanh(projection)
  modifier[!is.finite(projection)] <- 0
  projection_available <- is.finite(projection)
  gene_rna_support <- rc_gene_score(
    rna_logcpm, mode = "absolute", half_saturation = gene_half_saturation
  )
  gene_multiome_support <- .rc_integrate_regulatory_support(
    gene_rna_support, modifier, alpha = regulatory_alpha
  )
  reaction_expression <- rc_reaction_capacity(
    parsed, gene_multiome_support,
    promiscuity_mode = "none",
    and_method = gpr_and_method,
    or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  saturation <- do.call(rbind, lapply(seq_len(nrow(unit_meta)), function(i) {
    value <- projection[, i]
    finite <- is.finite(value)
    data.frame(
      metacell_id = units[[i]],
      condition = as.character(unit_meta[[condition_col]][[i]]),
      cell_type = as.character(unit_meta[[celltype_col]][[i]]),
      n_finite_targets = sum(finite),
      abs_projection_q50 = if (any(finite)) {
        stats::median(abs(value[finite]))
      } else {
        NA_real_
      },
      abs_projection_q90 = if (any(finite)) {
        unname(stats::quantile(abs(value[finite]), 0.9))
      } else {
        NA_real_
      },
      fraction_abs_gt_2 = if (any(finite)) mean(abs(value[finite]) > 2) else NA_real_,
      tanh_saturation_fraction = if (any(finite)) {
        mean(abs(tanh(value[finite])) > 0.95)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  list(
    schema_version = "regcompass_condition_grn_layer1_v3",
    reaction_expression = reaction_expression,
    rna_metacell_logcpm = rna_logcpm,
    gene_support_rna = gene_rna_support,
    gene_projection_absolute = projection_absolute,
    gene_projection_shared = projection_shared,
    gene_projection_deviation = projection_deviation,
    gene_projection_raw = projection,
    gene_regulatory_reliability = reliability,
    gene_regulatory_reliability_available = reliability_available,
    gene_regulatory_available = projection_available,
    gene_regulatory_modifier = modifier,
    gene_regulatory_model_projection = projection,
    gene_support_multiome = gene_multiome_support,
    projection_coverage = do.call(rbind, coverage),
    projection_saturation = saturation,
    parsed_gpr = parsed,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, rownames(rna_logcpm)),
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "SuperCell_condition_by_broad_cell_type",
    capacity_params = list(
      regulatory_alpha = regulatory_alpha,
      gene_half_saturation = gene_half_saturation,
      regulatory_mode = paste0(
        "Pando_cell_first_", projection_component, "_", comparison_support
      ),
      regulatory_reliability = "sqrt(max(0, pooled_sample_blocked_OOF_R2))",
      promiscuity_mode = "none",
      and_method = gpr_and_method,
      or_method = "sum",
      parallel = parallel
    ),
    projection_provenance = list(
      pando_schema = "pando_condition_grn_fit_v4",
      pando_version = grn_result$pando_installed_version,
      pando_file_fingerprint = grn_result$pando_file_fingerprint,
      supercell_membership = "misc$membership_table(cell_id, metacell_id)",
      projection_order = paste(
        "single-cell TF*ATAC -> pooled transform -> condition coefficient ->",
        "target sum -> SuperCell membership mean"
      ),
      comparison_support_requested = comparison_support,
      comparison_support_resolved = unique(
        unlist(lapply(coverage, function(x) x$comparison_support))
      ),
      component = projection_component
    )
  )
}

#' Project condition GRN effects and build RNA+ATAC reaction support
#'
#' Uses Pando's fitted single-cell TF-by-ATAC predictors and exact SuperCell2
#' cell-to-metacell membership. Interactions are never reconstructed from
#' metacell means.
#'
#' @export
rc_regcompass_step_layer1 <- function(
    grn, metacells, meta_modules, gem, outdir,
    projection_component = c("condition", "shared", "deviation"),
    comparison_support = c(
      "auto", "pairwise_common", "global_common",
      "condition_estimable", "strict"
    ),
    projection_mode = "metacell_specific",
    regulatory_reliability = "sample_blocked_oof",
    regulatory_alpha = 1,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer1", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  projection_component <- match.arg(projection_component)
  comparison_support <- match.arg(comparison_support)
  gpr_and_method <- match.arg(gpr_and_method)
  if (!identical(projection_mode, "metacell_specific") ||
      !identical(regulatory_reliability, "sample_blocked_oof")) {
    stop(
      "Canonical Layer 1 requires projection_mode='metacell_specific' and ",
      "regulatory_reliability='sample_blocked_oof'.", call. = FALSE
    )
  }
  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!identical(metacells$params, meta_modules$workflow_params)) {
    stop(
      "Metacell and meta-module stages use different workflow settings.",
      call. = FALSE
    )
  }
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  params <- metacells$params
  .rc_require_stage_gem(grn, gem, "grn")
  layer1 <- .rc_cell_first_projection_layer1(
    grn_result = grn$grn_result,
    metacell_object = metacells$metacell_object,
    membership = metacells$pooled$membership,
    gem = gem,
    metacell_meta = metacells$pooled$metacell_meta,
    sample_col = params$sample_col,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    projection_component = projection_component,
    comparison_support = comparison_support,
    regulatory_alpha = regulatory_alpha,
    gpr_and_method = gpr_and_method,
    gene_half_saturation = gene_half_saturation,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  layer1$workflow_params <- params
  layer1$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  class(layer1) <- c("regcompass_layer1_step", "list")
  .rc_validate_layer1_stage(
    layer1,
    workflow_params = params,
    gem = gem,
    argument = "layer1"
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer1 <- .rc_step_monitor_finish(layer1, monitor)
  saveRDS(layer1, file.path(outdir, "step_layer1.rds"))
  layer1
}
