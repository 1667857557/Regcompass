# Cell-type-specific routing between condition and standard Pando.

.rc_bind_pando_field <- function(results, field) {
  values <- lapply(results, function(x) x[[field]])
  values <- values[vapply(values, is.data.frame, logical(1))]
  values <- values[vapply(values, nrow, integer(1)) > 0L]
  if (!length(values)) return(data.frame())
  .rc_bind_frames_fill(values)
}

.rc_overlay_condition_fit_diagnostics <- function(
    all_edges, condition_fits, condition_col, celltype_col) {
  if (!nrow(all_edges) || !length(condition_fits)) return(all_edges)
  rows <- lapply(condition_fits, function(fit) {
    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    required <- c(
      "edge_id", "condition", "fit_status", "target_rsq", "target_rsq_oof",
      "target_model_evaluated", "target_model_supported",
      "target_rsq_threshold"
    )
    if (!all(required %in% colnames(coefficient))) {
      stop("Condition fit is missing RegCompass target-quality diagnostics.",
           call. = FALSE)
    }
    data.frame(
      edge_id = as.character(coefficient$edge_id),
      condition_value = as.character(coefficient$condition),
      celltype_value = as.character(fit$cell_type),
      fit_status = as.character(coefficient$fit_status),
      target_rsq = as.numeric(coefficient$target_rsq),
      target_rsq_oof = as.numeric(coefficient$target_rsq_oof),
      target_model_evaluated = as.logical(coefficient$target_model_evaluated),
      target_model_supported = as.logical(coefficient$target_model_supported),
      target_rsq_threshold = as.numeric(coefficient$target_rsq_threshold),
      stringsAsFactors = FALSE
    )
  })
  diagnostics <- do.call(rbind, rows)
  key <- paste(
    diagnostics$edge_id, diagnostics$condition_value,
    diagnostics$celltype_value, sep = "\001"
  )
  if (anyDuplicated(key)) {
    stop("Condition target-quality diagnostics contain duplicate edge keys.",
         call. = FALSE)
  }
  condition_rows <- rep(TRUE, nrow(all_edges))
  if ("analysis_mode" %in% colnames(all_edges)) {
    route <- as.character(all_edges$analysis_mode)
    condition_rows <- is.na(route) | !nzchar(route) | route == "condition_grn"
  }
  if (!any(condition_rows)) return(all_edges)
  required_edges <- c("edge_id", condition_col, celltype_col)
  if (!all(required_edges %in% colnames(all_edges))) {
    stop("Condition edge table lacks edge/condition/cell-type routing keys.",
         call. = FALSE)
  }
  edge_key <- paste(
    as.character(all_edges$edge_id[condition_rows]),
    as.character(all_edges[[condition_col]][condition_rows]),
    as.character(all_edges[[celltype_col]][condition_rows]), sep = "\001"
  )
  index <- match(edge_key, key)
  if (anyNA(index)) {
    stop("Condition edges cannot be aligned to target-quality diagnostics.",
         call. = FALSE)
  }
  fields <- c(
    "fit_status", "target_rsq", "target_rsq_oof",
    "target_model_evaluated", "target_model_supported", "target_rsq_threshold"
  )
  for (field in fields) {
    if (!field %in% colnames(all_edges)) all_edges[[field]] <- NA
    all_edges[[field]][condition_rows] <- diagnostics[[field]][index]
  }
  all_edges
}

.rc_merge_pando_results <- function(
    condition_result = NULL, standard_results = list(),
    condition_types = character(), standard_types = character(),
    condition_col, celltype_col, outdir,
    target_rsq_threshold = .RC_PANDO_TARGET_RSQ_THRESHOLD_DEFAULT) {
  target_rsq_threshold <- .rc_target_rsq_threshold(target_rsq_threshold)
  .rc_validate_pando_route_partition(condition_types, standard_types)
  results <- c(
    if (is.null(condition_result)) list() else list(condition_result),
    standard_results
  )
  if (!length(results)) stop("No Pando result was produced.", call. = FALSE)
  .rc_validate_pando_result_cell_partition(results)

  standard_objects <- list()
  for (value in standard_results) {
    for (name in names(value$standard_pando_objects)) {
      if (name %in% names(standard_objects)) {
        stop("Duplicated standard-Pando cell type: ", name, call. = FALSE)
      }
      standard_objects[[name]] <- value$standard_pando_objects[[name]]
    }
  }

  condition_fits <- if (is.null(condition_result)) {
    list()
  } else {
    condition_result$condition_grn_fits
  }
  if (length(condition_fits)) {
    condition_fits <- lapply(
      condition_fits,
      .rc_apply_condition_penalty_gate,
      target_rsq_threshold = target_rsq_threshold
    )
  }

  condition_objects <- if (is.null(condition_result)) {
    list()
  } else if (is.list(condition_result$pando_grn_data_by_cell_type) &&
             length(condition_result$pando_grn_data_by_cell_type)) {
    condition_result$pando_grn_data_by_cell_type
  } else if (inherits(condition_result$pando_grn_data, "GRNData") &&
             length(condition_fits) == 1L) {
    type <- as.character(condition_fits[[1L]]$cell_type)[[1L]]
    stats::setNames(list(condition_result$pando_grn_data), type)
  } else {
    list()
  }
  if (length(condition_fits) &&
      !setequal(
        vapply(condition_fits, function(fit) as.character(fit$cell_type)[[1L]],
               character(1)),
        names(condition_objects)
      )) {
    stop("Condition Pando objects and fit contracts cover different cell types.",
         call. = FALSE)
  }

  paired_meta <- .rc_bind_pando_field(results, "paired_cell_metadata")
  if (nrow(paired_meta)) {
    paired_meta <- paired_meta[!duplicated(paired_meta$cell_id), , drop = FALSE]
  }
  routing <- rbind(
    if (length(condition_types)) data.frame(
      cell_type = condition_types, analysis_mode = "condition_grn",
      stringsAsFactors = FALSE
    ),
    if (length(standard_types)) data.frame(
      cell_type = standard_types, analysis_mode = "standard_pando",
      stringsAsFactors = FALSE
    )
  )
  rownames(routing) <- NULL
  mode <- if (length(condition_types) && length(standard_types)) {
    "mixed_pando"
  } else if (length(condition_types)) {
    "condition_grn"
  } else {
    "standard_pando"
  }

  all_edges <- .rc_bind_pando_field(results, "tf_peak_gene_condition_all")
  all_edges <- .rc_overlay_condition_fit_diagnostics(
    all_edges, condition_fits, condition_col, celltype_col
  )
  if (nrow(all_edges)) {
    route <- if ("analysis_mode" %in% colnames(all_edges)) {
      as.character(all_edges$analysis_mode)
    } else {
      rep(if (length(condition_fits)) "condition_grn" else "standard_pando",
          nrow(all_edges))
    }
    condition_rows <- is.na(route) | !nzchar(route) | route == "condition_grn"
    standard_rows <- route == "standard_pando"
    all_edges$penalty_eligible <- FALSE
    all_edges$active_in_condition <- FALSE

    if (any(condition_rows)) {
      gate <- .rc_condition_penalty_gate(
        all_edges[condition_rows, , drop = FALSE],
        target_rsq_threshold = target_rsq_threshold
      )
      all_edges$penalty_eligible[condition_rows] <- gate
      all_edges$active_in_condition[condition_rows] <- gate
    }
    if (any(standard_rows)) {
      annotated <- .rc_annotate_standard_pando_edges(
        all_edges[standard_rows, , drop = FALSE],
        target_rsq_threshold = target_rsq_threshold
      )
      for (field in c(
        "target_rsq", "target_rsq_oof", "target_model_evaluated",
        "target_model_supported", "target_rsq_threshold"
      )) {
        if (!field %in% colnames(all_edges)) all_edges[[field]] <- NA
        all_edges[[field]][standard_rows] <- annotated[[field]]
      }
      estimate <- suppressWarnings(as.numeric(annotated$estimate))
      padj <- suppressWarnings(as.numeric(annotated$padj))
      standard_gate <- is.finite(estimate) & is.finite(padj) &
        padj < as.numeric(annotated$padj_threshold %||% 0.05) &
        annotated$target_model_supported %in% TRUE
      if ("estimable" %in% colnames(annotated)) {
        standard_gate <- standard_gate & annotated$estimable %in% TRUE
      }
      all_edges$penalty_eligible[standard_rows] <- standard_gate
      all_edges$active_in_condition[standard_rows] <- standard_gate
    }
  }
  active_edges <- if (nrow(all_edges)) {
    all_edges[all_edges$penalty_eligible %in% TRUE, , drop = FALSE]
  } else {
    data.frame()
  }

  answer <- list(
    schema_version = "regcompass_celltype_routed_pando",
    analysis_mode = mode,
    cell_type_analysis_mode = routing,
    condition_coefficients_calculated = length(condition_fits) > 0L,
    pando_fit_schema = if (length(condition_fits)) {
      .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
    } else NA_character_,
    pando_model_schema = if (length(condition_fits)) {
      .RC_PANDO_CONDITION_GRN_MODEL_SCHEMA
    } else NA_character_,
    pando_installed_version = as.character(utils::packageVersion("Pando")),
    pando_grn_data = if (length(condition_objects) == 1L) {
      condition_objects[[1L]]
    } else {
      NULL
    },
    pando_grn_data_by_cell_type = condition_objects,
    pando_object_scope = list(
      cell_types = names(condition_objects),
      preserves_cell_type_peak_space = TRUE,
      combined_grndata = length(condition_objects) == 1L
    ),
    standard_pando_objects = standard_objects,
    condition_grn_fits = condition_fits,
    target_metabolic_genes = unique(unlist(
      lapply(results, `[[`, "target_metabolic_genes"), use.names = FALSE
    )),
    condition_fit_status = .rc_bind_pando_field(results, "condition_fit_status"),
    pando_network_index = if (is.null(condition_result)) data.frame() else
      condition_result$pando_network_index,
    pando_fit_diagnostics = if (is.null(condition_result)) data.frame() else
      condition_result$pando_fit_diagnostics,
    pando_execution_summary = if (is.null(condition_result)) list() else
      condition_result$pando_execution_summary,
    tf_peak_gene_universal = .rc_bind_pando_field(
      results, "tf_peak_gene_universal"
    ),
    tf_peak_gene_condition_all = all_edges,
    tf_peak_gene_condition = active_edges,
    tf_peak_gene_condition_effect_all = if (is.null(condition_result)) {
      data.frame()
    } else {
      all_edges[
        if ("analysis_mode" %in% colnames(all_edges)) {
          is.na(all_edges$analysis_mode) |
            !nzchar(as.character(all_edges$analysis_mode)) |
            as.character(all_edges$analysis_mode) == "condition_grn"
        } else rep(TRUE, nrow(all_edges)),
        , drop = FALSE
      ]
    },
    tf_peak_gene_condition_effect = if (is.null(condition_result)) {
      data.frame()
    } else if (nrow(active_edges) &&
               "analysis_mode" %in% colnames(active_edges)) {
      active_edges[
        is.na(active_edges$analysis_mode) |
          !nzchar(as.character(active_edges$analysis_mode)) |
          as.character(active_edges$analysis_mode) == "condition_grn",
        , drop = FALSE
      ]
    } else {
      active_edges
    },
    tf_peak_gene_condition_contrasts = .rc_bind_pando_field(
      results, "tf_peak_gene_condition_contrasts"
    ),
    paired_cell_ids = unique(as.character(paired_meta$cell_id)),
    paired_cell_metadata = paired_meta,
    target_rsq_threshold = target_rsq_threshold,
    target_model_quality_policy =
      "selected-lambda final full-data R2 gate; OOF R2 diagnostic only",
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF",
      cell_type_routing = paste(
        "condition GRN for at least two retained conditions;",
        "standard Pando otherwise"
      ),
      condition_effect_filter = paste(
        "consume Pando condition-specific estimable BH-active ridge edge;",
        "do not redefine Pando significance; then require target fit_status ok",
        "and selected-lambda final full-data R2 >= target_rsq_threshold"
      ),
      standard_edge_filter = paste(
        "estimable when available and adjusted P below configured threshold;",
        "require selected-lambda final full-data R2 >= target_rsq_threshold"
      ),
      projection =
        "beta times metacell-mean TF times metacell-mean ATAC"
    ),
    group_cols = c(condition_col, celltype_col)
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(
    answer$condition_fit_status,
    file.path(outdir, "pando_group_status.tsv.gz")
  )
  .rc_write_tsv_gz(
    answer$tf_peak_gene_condition_all,
    file.path(outdir, "pando_tf_peak_gene_all.tsv.gz")
  )
  .rc_write_tsv_gz(
    answer$tf_peak_gene_condition,
    file.path(outdir, "pando_tf_peak_gene_active.tsv.gz")
  )
  if (nrow(answer$tf_peak_gene_condition_contrasts)) {
    .rc_write_tsv_gz(
      answer$tf_peak_gene_condition_contrasts,
      file.path(outdir, "pando_tf_peak_gene_condition_contrasts.tsv.gz")
    )
  }
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_active_target_penalty_q <- function(
    grn_result, unit_meta, condition_col, celltype_col, template) {
  q <- as.matrix(template)
  q[,] <- NA_real_
  all_edges <- grn_result$tf_peak_gene_condition_all
  if (!is.data.frame(all_edges) || !nrow(all_edges)) return(q)
  required <- c(
    "target", condition_col, celltype_col,
    "target_model_evaluated", "penalty_eligible"
  )
  if (!all(required %in% colnames(all_edges)) ||
      !all(c("unit_id", condition_col, celltype_col) %in% colnames(unit_meta))) {
    stop(
      "Target penalty q requires evaluated-target and penalty eligibility ",
      "labels aligned to metacell metadata.", call. = FALSE
    )
  }
  target <- tolower(trimws(as.character(all_edges$target)))
  condition <- trimws(as.character(all_edges[[condition_col]]))
  cell_type <- trimws(as.character(all_edges[[celltype_col]]))
  if (anyNA(target) || any(!nzchar(target)) ||
      anyNA(condition) || any(!nzchar(condition)) ||
      anyNA(cell_type) || any(!nzchar(cell_type))) {
    stop("Pando edges contain incomplete penalty-routing labels.",
         call. = FALSE)
  }
  all_edges$target <- target
  all_edges[[condition_col]] <- condition
  all_edges[[celltype_col]] <- cell_type
  strata <- unique(all_edges[, c(condition_col, celltype_col), drop = FALSE])
  for (i in seq_len(nrow(strata))) {
    condition_value <- as.character(strata[[condition_col]][[i]])
    celltype_value <- as.character(strata[[celltype_col]][[i]])
    edge_keep <- all_edges[[condition_col]] == condition_value &
      all_edges[[celltype_col]] == celltype_value
    units <- intersect(
      as.character(unit_meta$unit_id[
        as.character(unit_meta[[condition_col]]) == condition_value &
          as.character(unit_meta[[celltype_col]]) == celltype_value
      ]),
      colnames(q)
    )
    if (!length(units)) next
    evaluated_targets <- intersect(
      unique(as.character(all_edges$target[
        edge_keep & all_edges$target_model_evaluated %in% TRUE
      ])),
      rownames(q)
    )
    if (length(evaluated_targets)) {
      q[evaluated_targets, units] <- 0
    }
    active_targets <- intersect(
      unique(as.character(all_edges$target[
        edge_keep & all_edges$penalty_eligible %in% TRUE
      ])),
      rownames(q)
    )
    if (length(active_targets)) {
      q[active_targets, units] <- 1
    }
  }
  q
}

.rc_project_pando_by_celltype <- function(
    grn_result, membership, unit_meta, genes,
    condition_col, celltype_col, rna_assay, atac_assay) {
  template <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  projection <- template
  coverage <- list()
  origins <- schemas <- projection_names <- policies <- character()
  if (length(grn_result$condition_grn_fits)) {
    part <- .rc_condition_pando_projection(
      grn_result, membership, unit_meta, genes, rna_assay, atac_assay
    )
    projection <- .rc_overlay_projection(
      projection, part$projection, "Condition-Pando"
    )
    coverage[[length(coverage) + 1L]] <- part$coverage
    origins <- c(origins, part$origin)
    schemas <- c(schemas, part$pando_schema)
    projection_names <- c(projection_names, part$projection_name)
    policies <- c(policies, part$nonestimable_policy)
  }
  if (length(grn_result$standard_pando_objects)) {
    standard <- .rc_standard_pando_projection(
      grn_result, membership, unit_meta, condition_col, celltype_col,
      rna_assay = rna_assay, atac_assay = atac_assay,
      target_genes = genes
    )
    projection <- .rc_overlay_projection(
      projection, standard$projection, "Standard-Pando"
    )
    coverage[[length(coverage) + 1L]] <- standard$coverage
    origins <- c(origins, standard$projection_origin)
    schemas <- c(schemas, "standard_pando_network")
    projection_names <- c(projection_names, "standard_pando_full_fit")
    policies <- c(policies, "not_applicable_standard_pando")
  }
  if (!length(origins)) stop("No Pando projection route is available.", call. = FALSE)

  reliability <- .rc_active_target_penalty_q(
    grn_result = grn_result,
    unit_meta = unit_meta,
    condition_col = condition_col,
    celltype_col = celltype_col,
    template = template
  )
  coverage_table <- .rc_bind_frames_fill(coverage)
  if (nrow(coverage_table)) {
    coverage_table$penalty_q_definition <- paste(
      "q=1 for target with penalty-eligible edge; q=0 for evaluated target",
      "without eligible edge; q=NA for unavailable/not-evaluated target;",
      "final full-data R2 is a target gate and OOF R2 is diagnostic only"
    )
  }
  list(
    projection = projection,
    reliability = reliability,
    coverage = coverage_table,
    origin = paste(unique(origins), collapse = ";"),
    pando_schema = paste(unique(schemas), collapse = ";"),
    projection_name = paste(unique(projection_names), collapse = ";"),
    nonestimable_policy = paste(unique(policies), collapse = ";"),
    condition_coefficients_calculated =
      length(grn_result$condition_grn_fits) > 0L,
    cell_type_analysis_mode = grn_result$cell_type_analysis_mode,
    penalty_q_definition =
      "q=1 eligible; q=0 evaluated neutral; q=NA unavailable; OOF R2 diagnostic"
  )
}
