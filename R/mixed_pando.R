# Cell-type-specific routing between common-dictionary and standard Pando.

.rc_bind_pando_field <- function(results, field) {
  values <- lapply(results, function(x) x[[field]])
  values <- values[vapply(values, is.data.frame, logical(1))]
  values <- values[vapply(values, nrow, integer(1)) > 0L]
  if (!length(values)) return(data.frame())
  .rc_bind_frames_fill(values)
}

.rc_merge_pando_results <- function(
    condition_result = NULL, standard_results = list(),
    condition_types = character(), standard_types = character(),
    condition_col, celltype_col, outdir) {
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
    condition_fits <- lapply(condition_fits, .rc_apply_condition_penalty_gate)
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
  active_edges <- .rc_bind_pando_field(results, "tf_peak_gene_condition")
  if (nrow(all_edges) && length(condition_fits)) {
    condition_rows <- rep(TRUE, nrow(all_edges))
    if ("analysis_mode" %in% colnames(all_edges)) {
      route <- as.character(all_edges$analysis_mode)
      condition_rows <- is.na(route) | !nzchar(route) | route == "condition_grn"
    }
    if (any(condition_rows)) {
      gate <- .rc_condition_penalty_gate(
        all_edges[condition_rows, , drop = FALSE]
      )
      all_edges$significant[condition_rows] <- gate
      all_edges$penalty_effect[condition_rows] <- ifelse(
        gate,
        suppressWarnings(as.numeric(all_edges$estimate[condition_rows])),
        0
      )
      if ("penalty_eligible" %in% colnames(all_edges)) {
        all_edges$penalty_eligible[condition_rows] <- gate
      }
      if ("active_in_condition" %in% colnames(all_edges)) {
        all_edges$active_in_condition[condition_rows] <- gate
      }
      condition_active <- all_edges[condition_rows &
        all_edges$significant %in% TRUE, , drop = FALSE]
      standard_active <- if (nrow(active_edges) &&
          "analysis_mode" %in% colnames(active_edges)) {
        active_edges[
          as.character(active_edges$analysis_mode) == "standard_pando",
          , drop = FALSE
        ]
      } else {
        data.frame()
      }
      active_edges <- .rc_bind_frames_fill(list(
        condition_active, standard_active
      ))
    }
  }

  answer <- list(
    schema_version = "regcompass_celltype_routed_pando",
    analysis_mode = mode,
    cell_type_analysis_mode = routing,
    condition_coefficients_calculated = length(condition_fits) > 0L,
    pando_fit_schema = if (length(condition_fits)) {
      .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
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
    paired_cell_ids = unique(as.character(paired_meta$cell_id)),
    paired_cell_metadata = paired_meta,
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF",
      cell_type_routing = paste(
        "condition GRN for at least two retained conditions;",
        "standard Pando otherwise"
      ),
      condition_effect_filter = "estimable and BH adjusted P below 0.05",
      standard_edge_filter =
        "estimable when available and adjusted P below 0.05",
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
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_project_pando_by_celltype <- function(
    grn_result, membership, unit_meta, genes,
    condition_col, celltype_col, rna_assay, atac_assay) {
  template <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  projection <- reliability <- template
  coverage <- list()
  origins <- schemas <- projection_names <- policies <- character()
  if (length(grn_result$condition_grn_fits)) {
    part <- .rc_condition_pando_projection(
      grn_result, membership, unit_meta, genes, rna_assay, atac_assay
    )
    projection <- .rc_overlay_projection(
      projection, part$projection, "Condition-Pando"
    )
    reliability <- .rc_overlay_projection(
      reliability, part$reliability, "Condition-Pando reliability"
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
    reliability <- .rc_overlay_projection(
      reliability, standard$reliability, "Standard-Pando reliability"
    )
    coverage[[length(coverage) + 1L]] <- standard$coverage
    origins <- c(origins, standard$projection_origin)
    schemas <- c(schemas, "standard_pando_network")
    projection_names <- c(projection_names, "standard_pando_full_fit")
    policies <- c(policies, "not_applicable_standard_pando")
  }
  if (!length(origins)) stop("No Pando projection route is available.", call. = FALSE)
  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = paste(unique(origins), collapse = ";"),
    pando_schema = paste(unique(schemas), collapse = ";"),
    projection_name = paste(unique(projection_names), collapse = ";"),
    nonestimable_policy = paste(unique(policies), collapse = ";"),
    condition_coefficients_calculated =
      length(grn_result$condition_grn_fits) > 0L,
    cell_type_analysis_mode = grn_result$cell_type_analysis_mode
  )
}
