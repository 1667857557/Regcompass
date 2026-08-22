# Cell-type-specific routing between condition and standard Pando.

.rc_bind_pando_field <- function(results, field) {
  values <- lapply(results, function(x) x[[field]])
  values <- values[vapply(values, is.data.frame, logical(1))]
  values <- values[vapply(values, nrow, integer(1)) > 0L]
  if (!length(values)) return(data.frame())
  .rc_bind_frames_fill(values)
}

.rc_condition_penalty_gate_by_celltype <- function(coefficient, celltype_col) {
  if (!is.data.frame(coefficient) || !nrow(coefficient)) return(logical())
  if (!is.character(celltype_col) || length(celltype_col) != 1L ||
      is.na(celltype_col) || !nzchar(trimws(celltype_col)) ||
      !celltype_col %in% colnames(coefficient)) {
    stop("Condition penalty gating requires its broad cell-type column.",
         call. = FALSE)
  }
  cell_type <- trimws(as.character(coefficient[[celltype_col]]))
  if (anyNA(cell_type) || any(!nzchar(cell_type))) {
    stop("Condition penalty gating received incomplete broad cell types.",
         call. = FALSE)
  }
  rows <- unname(split(
    seq_len(nrow(coefficient)), factor(cell_type, levels = unique(cell_type))
  ))
  gate <- logical(nrow(coefficient))
  for (index in rows) {
    gate[index] <- .rc_condition_penalty_gate(
      coefficient[index, , drop = FALSE]
    )
  }
  gate
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
      # E-star owns continuous fixed-dictionary production effects. Pando owns
      # the common exact-edge whole-network BH topology. RegCompass validates
      # that topology and never rebuilds a condition-local edge-selection rule.
      condition_table <- all_edges[condition_rows, , drop = FALSE]
      gate <- .rc_condition_penalty_gate_by_celltype(
        condition_table, celltype_col = celltype_col
      )
      target_rsq <- .rc_condition_target_rsq(condition_table)
      rsq_threshold <- .rc_target_rsq_threshold()
      if (!"target_rsq" %in% colnames(all_edges)) {
        all_edges$target_rsq <- suppressWarnings(as.numeric(all_edges$rsq))
      }
      all_edges$target_rsq[condition_rows] <- target_rsq
      if (!"target_rsq_threshold" %in% colnames(all_edges)) {
        all_edges$target_rsq_threshold <- rsq_threshold
      }
      all_edges$target_rsq_threshold[condition_rows] <- rsq_threshold
      if (!"target_model_supported" %in% colnames(all_edges)) {
        all_edges$target_model_supported <- FALSE
      }
      fit_status <- trimws(as.character(condition_table$fit_status))
      all_edges$target_model_supported[condition_rows] <-
        !is.na(fit_status) & fit_status == "ok" &
        is.finite(target_rsq) & target_rsq >= rsq_threshold
      if (!"penalty_eligible" %in% colnames(all_edges)) {
        all_edges$penalty_eligible <- FALSE
      }
      all_edges$penalty_eligible[condition_rows] <- gate
      if (!"active_in_condition" %in% colnames(all_edges)) {
        all_edges$active_in_condition <- FALSE
      }
      all_edges$active_in_condition[condition_rows] <- gate
      condition_active <- all_edges[
        condition_rows & all_edges$penalty_eligible %in% TRUE,
        , drop = FALSE
      ]
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
    schema_version = "regcompass_celltype_routed_pando_v2",
    analysis_mode = mode,
    cell_type_analysis_mode = routing,
    condition_coefficients_calculated = length(condition_fits) > 0L,
    target_rsq_threshold = .rc_target_rsq_threshold(),
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
    pando_edge_inference = .rc_bind_pando_field(results, "pando_edge_inference"),
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
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF",
      cell_type_routing = paste(
        "condition GRN for at least two retained conditions;",
        "standard Pando otherwise"
      ),
      condition_effect_filter = paste(
        "consume Pando exact-edge whole-network BH common topology and retain",
        "each condition's finite continuous E-star z=0.25 beta_E; target R2",
        "is diagnostic only"
      ),
      condition_inference = paste(
        "condition-local no-fusion Gaussian coefficient tests -> one exact-edge",
        "omnibus P -> whole-cell-type-network BH"
      ),
      condition_contrast_filter =
        "biological differential claims require contrast_identifiable == TRUE",
      standard_edge_filter = paste(
        "estimable when available, adjusted P below the configured threshold,",
        paste0("and selected-lambda full-data R2 >= ",
               format(.rc_target_rsq_threshold(), trim = TRUE))
      ),
      projection =
        "continuous beta_E times canonical RegCompass metacell TF-ATAC exposure"
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
  if (nrow(answer$pando_edge_inference)) {
    .rc_write_tsv_gz(
      answer$pando_edge_inference,
      file.path(outdir, "pando_exact_edge_inference.tsv.gz")
    )
  }
  if (nrow(answer$tf_peak_gene_condition_contrasts)) {
    .rc_write_tsv_gz(
      answer$tf_peak_gene_condition_contrasts,
      file.path(outdir, "pando_tf_peak_gene_condition_contrasts.tsv.gz")
    )
  }
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_validate_penalty_q_table <- function(
    table, unit_meta, condition_col, celltype_col, label) {
  if (!is.data.frame(table) || !nrow(table)) return(table)
  required <- c("target", condition_col, celltype_col)
  if (!all(required %in% colnames(table)) ||
      !all(c("unit_id", condition_col, celltype_col) %in% colnames(unit_meta))) {
    stop(
      label, " penalty routing requires target, condition and cell-type labels ",
      "aligned to metacell metadata.", call. = FALSE
    )
  }
  table$target <- tolower(trimws(as.character(table$target)))
  table[[condition_col]] <- trimws(as.character(table[[condition_col]]))
  table[[celltype_col]] <- trimws(as.character(table[[celltype_col]]))
  if (anyNA(table$target) || any(!nzchar(table$target)) ||
      anyNA(table[[condition_col]]) || any(!nzchar(table[[condition_col]])) ||
      anyNA(table[[celltype_col]]) || any(!nzchar(table[[celltype_col]]))) {
    stop(label, " Pando edges contain incomplete penalty-routing labels.",
         call. = FALSE)
  }
  table
}

.rc_penalty_evaluated_rows <- function(table) {
  if (!is.data.frame(table) || !nrow(table)) return(logical())
  condition_route <- if ("analysis_mode" %in% colnames(table)) {
    route <- as.character(table$analysis_mode)
    is.na(route) | !nzchar(route) | route == "condition_grn"
  } else {
    rep(FALSE, nrow(table))
  }
  fit_ok <- if ("fit_status" %in% colnames(table)) {
    status <- trimws(as.character(table$fit_status))
    !is.na(status) & status == "ok"
  } else {
    rep(TRUE, nrow(table))
  }
  # Conditional targets are evaluated from production fit status and finite
  # beta_E. Target R2 is not an eligibility gate. Standard Pando retains the
  # legacy finite-R2 evaluation semantics.
  rsq <- if ("target_rsq" %in% colnames(table)) {
    suppressWarnings(as.numeric(table$target_rsq))
  } else if ("rsq" %in% colnames(table)) {
    suppressWarnings(as.numeric(table$rsq))
  } else {
    rep(NA_real_, nrow(table))
  }
  standard_evaluated <- fit_ok & is.finite(rsq)
  condition_evaluated <- fit_ok
  if ("penalty_effect" %in% colnames(table)) {
    condition_evaluated <- condition_evaluated &
      is.finite(suppressWarnings(as.numeric(table$penalty_effect)))
  }
  ifelse(condition_route, condition_evaluated, standard_evaluated)
}

.rc_set_penalty_q_by_stratum <- function(
    q, table, unit_meta, condition_col, celltype_col, value) {
  if (!nrow(table)) return(q)
  strata <- unique(table[, c(condition_col, celltype_col), drop = FALSE])
  for (i in seq_len(nrow(strata))) {
    condition_value <- as.character(strata[[condition_col]][[i]])
    celltype_value <- as.character(strata[[celltype_col]][[i]])
    edge_keep <- table[[condition_col]] == condition_value &
      table[[celltype_col]] == celltype_value
    targets <- intersect(unique(as.character(table$target[edge_keep])), rownames(q))
    units <- intersect(
      as.character(unit_meta$unit_id[
        as.character(unit_meta[[condition_col]]) == condition_value &
          as.character(unit_meta[[celltype_col]]) == celltype_value
      ]),
      colnames(q)
    )
    if (length(targets) && length(units)) q[targets, units] <- value
  }
  q
}

.rc_active_target_penalty_q <- function(
    grn_result, unit_meta, condition_col, celltype_col, template) {
  q <- as.matrix(template)
  q[,] <- NA_real_
  evaluated <- .rc_validate_penalty_q_table(
    grn_result$tf_peak_gene_condition_all,
    unit_meta, condition_col, celltype_col, "Evaluated-target"
  )
  if (nrow(evaluated)) {
    evaluated <- evaluated[.rc_penalty_evaluated_rows(evaluated), , drop = FALSE]
  }
  active <- .rc_validate_penalty_q_table(
    grn_result$tf_peak_gene_condition,
    unit_meta, condition_col, celltype_col, "Active-target"
  )

  # Conditional route: q=0 for a validly fitted target with no admitted common
  # exact edge; q=1 when at least one whole-network-BH-supported exact edge is
  # available for projection; q=NA when the target was not validly evaluated.
  # Target R2 remains diagnostic only. Standard Pando retains its existing
  # filtered-edge three-state semantics.
  q <- .rc_set_penalty_q_by_stratum(
    q, evaluated, unit_meta, condition_col, celltype_col, 0
  )
  q <- .rc_set_penalty_q_by_stratum(
    q, active, unit_meta, condition_col, celltype_col, 1
  )
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
      "condition route: q=1 when a valid target has at least one common",
      "whole-network-BH-supported exact edge, q=0 when valid but unsupported,",
      "q=NA when unavailable; target R2 diagnostic only; standard Pando keeps",
      "its existing filtered-edge q semantics"
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
    penalty_q_definition = paste(
      "condition route uses the common exact-edge whole-network BH topology;",
      "target R2 is diagnostic only; standard Pando retains legacy filtering"
    )
  )
}

