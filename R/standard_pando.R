# Standard Pando fallback used when no multi-level condition is available.

.rc_standard_pando_padj_fixed <- 0.05
.rc_standard_pando_min_abs_fixed <- 0.01

.rc_standard_pando_cell_types <- function(metadata, celltype_col, cell_type) {
  .rc_validate_celltype_metadata(metadata, celltype_col)
  available <- unique(as.character(metadata[[celltype_col]]))
  selected <- if (is.null(cell_type)) available else unique(as.character(cell_type))
  missing <- setdiff(selected, available)
  if (length(missing)) {
    stop("Requested cell types were not found: ", paste(missing, collapse = ", "),
         ".", call. = FALSE)
  }
  selected
}

.rc_standard_pando_infer_args <- function(args) {
  if (!is.list(args)) stop("`pando_infer_args` must be a list.", call. = FALSE)
  forbidden <- intersect(names(args), c(
    "object", "genes", "network_name", "aggregate_rna_col",
    "aggregate_peaks_col", "parallel"
  ))
  if (length(forbidden)) {
    stop("Standard Pando inference arguments cannot override managed fields: ",
         paste(forbidden, collapse = ", "), ".", call. = FALSE)
  }
  condition_only <- intersect(names(args), c(
    "candidate_screen", "condition_mix", "condition_weight",
    "reference_condition", "comparison_conditions", "nlambda", "lambda",
    "lambda_min_ratio", "outer_nfolds", "inner_nfolds",
    "lambda_selection", "active_tol", "seed", "max_iter",
    "tol_objective", "tol_coef", "BPPARAM"
  ))
  requested_scale <- args$scale %||% NULL
  if (length(condition_only)) args[condition_only] <- NULL
  if (!is.null(args$interaction_term) &&
      !identical(as.character(args$interaction_term), ":")) {
    stop("Standard RegCompass projection requires `interaction_term = ':'`.",
         call. = FALSE)
  }
  args$scale <- FALSE
  answer <- modifyList(list(interaction_term = ":", scale = FALSE), args)
  attr(answer, "standard_fallback_adjustments") <- list(
    dropped_condition_arguments = condition_only,
    scale_forced_false = !identical(requested_scale, FALSE),
    reason = "one_effective_condition_uses_original_pando_projection_scale"
  )
  answer
}

.rc_filter_standard_pando_edges <- function(table) {
  required <- c("estimate", "padj")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop("Standard Pando requires estimate and padj columns.", call. = FALSE)
  }
  estimate <- suppressWarnings(as.numeric(table$estimate))
  padj <- suppressWarnings(as.numeric(table$padj))
  keep <- is.finite(estimate) &
    abs(estimate) > .rc_standard_pando_min_abs_fixed &
    is.finite(padj) & padj < .rc_standard_pando_padj_fixed
  answer <- table[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    padj = "< 0.05",
    absolute_estimate = paste0(
      "> ", format(.rc_standard_pando_min_abs_fixed, scientific = FALSE)
    )
  )
  answer
}

.rc_fit_standard_pando_by_cell_type <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col, celltype_col = "cell_type", cell_type = NULL,
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(), pando_infer_args = list(),
    save_pando_objects = TRUE, parallel = FALSE,
    progress_monitor = NULL,
    species = c("auto", "human", "mouse")) {
  .rc_step_monitor_event(
    progress_monitor, "standard_runtime_check",
    "validating standard Pando runtime and inputs", current = 5L
  )
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Package 'Pando' is required.", call. = FALSE)
  }
  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  .rc_validate_condition_celltype_metadata(
    object@meta.data, condition_col, celltype_col
  )
  cell_types <- .rc_standard_pando_cell_types(
    object@meta.data, celltype_col, cell_type
  )
  .rc_step_monitor_event(
    progress_monitor, "standard_design",
    "resolved standard Pando cell-type jobs", current = 5L,
    context = list(cell_types = length(cell_types))
  )
  condition_levels <- unique(as.character(object@meta.data[[condition_col]]))
  if (length(condition_levels) != 1L) {
    stop("Standard Pando mode requires exactly one effective condition level.",
         call. = FALSE)
  }
  infer_args <- .rc_standard_pando_infer_args(pando_infer_args)
  fallback_adjustments <- attr(
    infer_args, "standard_fallback_adjustments", exact = TRUE
  )
  attr(infer_args, "standard_fallback_adjustments") <- NULL
  if (length(fallback_adjustments$dropped_condition_arguments) ||
      isTRUE(fallback_adjustments$scale_forced_false)) {
    message(
      "One effective condition detected: using standard Pando with scale = FALSE",
      if (length(fallback_adjustments$dropped_condition_arguments)) {
        paste0(
          "; ignored condition-only arguments: ",
          paste(fallback_adjustments$dropped_condition_arguments, collapse = ", ")
        )
      } else {
        ""
      }
    )
  }
  if (is.null(pfm)) pfm <- .rc_default_pando_motifs()
  if (!"regions" %in% names(pando_initiate_args) ||
      is.null(pando_initiate_args$regions)) {
    pando_initiate_args$regions <- .rc_default_pando_regions(species)
  }
  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  target_upper <- intersect(toupper(rna_genes), toupper(metabolic_genes))
  target_genes <- rna_genes[toupper(rna_genes) %in% target_upper]
  if (!length(target_genes)) {
    stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)
  }
  .rc_step_monitor_event(
    progress_monitor, "standard_target_selection",
    "resolved metabolic target genes", current = 6L,
    context = list(targets = length(target_genes))
  )
  filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Standard Pando fallback"
  )
  object <- filtered$object
  .rc_step_monitor_event(
    progress_monitor, "standard_atac_filter",
    "removed globally zero ATAC features", current = 6L,
    context = list(
      cells = ncol(object),
      removed_atac_features = filtered$n_removed %||% NA_integer_
    )
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "pando_objects"), recursive = TRUE,
             showWarnings = FALSE)
  objects <- list()
  all_rows <- list()
  active_rows <- list()
  status_rows <- list()
  selected_cells <- character()
  for (cell_index in seq_along(cell_types)) {
    value <- cell_types[[cell_index]]
    job_context <- list(
      cell_type = value, job = cell_index, jobs = length(cell_types)
    )
    .rc_step_monitor_event(
      progress_monitor, "standard_cell_type_start",
      "starting standard Pando cell-type job", current = 7L,
      context = job_context
    )
    cells <- rownames(object@meta.data)[
      as.character(object@meta.data[[celltype_col]]) == value
    ]
    if (length(cells) < as.integer(min_cells)) {
      stop("Cell type `", value, "` contains fewer than ", min_cells,
           " cells.", call. = FALSE)
    }
    one <- subset(object, cells = cells)
    init <- list(object = one, peak_assay = atac_assay, rna_assay = rna_assay)
    init[names(pando_initiate_args)] <- NULL
    .rc_step_monitor_event(
      progress_monitor, "standard_candidate_initialization",
      "initializing regulatory candidate space", current = 7L,
      context = c(job_context, list(cells = length(cells)))
    )
    grn <- do.call(Pando::initiate_grn, c(init, pando_initiate_args))
    .rc_step_monitor_event(
      progress_monitor, "standard_candidate_initialization_complete",
      "regulatory candidate space initialized", current = 7L,
      context = c(job_context, list(cells = length(cells)))
    )
    motif_args <- .rc_regcompass_motif_args(pando_motif_args)
    motif <- list(object = grn, pfm = pfm, genome = genome)
    motif[names(motif_args)] <- NULL
    .rc_step_monitor_event(
      progress_monitor, "standard_motif_mapping",
      "mapping binary peak-by-motif candidates", current = 8L,
      context = job_context
    )
    grn <- do.call(Pando::find_motifs, c(motif, motif_args))
    .rc_step_monitor_event(
      progress_monitor, "standard_motif_mapping_complete",
      "completed binary peak-by-motif mapping", current = 8L,
      context = job_context
    )
    infer <- list(
      object = grn,
      genes = target_genes,
      network_name = "regcompass_standard_grn",
      parallel = isTRUE(parallel)
    )
    infer[names(infer_args)] <- NULL
    .rc_step_monitor_event(
      progress_monitor, "standard_grn_fit",
      "fitting standard Pando GRN", current = 9L,
      context = c(job_context, list(targets = length(target_genes)))
    )
    grn <- do.call(Pando::infer_grn, c(infer, infer_args))
    .rc_step_monitor_event(
      progress_monitor, "standard_grn_fit_complete",
      "standard Pando GRN fit completed", current = 9L,
      context = job_context
    )
    .rc_step_monitor_event(
      progress_monitor, "standard_contract_extraction",
      "extracting standard Pando edge contract", current = 10L,
      context = job_context
    )
    extracted <- rc_extract_pando_tf_peak_gene(
      grn_object = grn,
      sample_id = value,
      padj_threshold = .rc_standard_pando_padj_fixed,
      require_padj = TRUE
    )
    extracted$significant <- .rc_filter_standard_pando_edges(extracted$all)
    .rc_step_monitor_event(
      progress_monitor, "standard_contract_extraction_complete",
      "extracted standard Pando edge contract", current = 10L,
      context = c(job_context, list(
        all_edges = nrow(extracted$all),
        active_edges = nrow(extracted$significant)
      ))
    )
    add_design <- function(tab) {
      if (!nrow(tab)) return(tab)
      tab[[condition_col]] <- condition_levels[[1L]]
      tab[[celltype_col]] <- value
      tab$group_id <- rc_make_stratum_id(
        tab[1L, c(condition_col, celltype_col), drop = FALSE],
        c(condition_col, celltype_col)
      )
      tab$effect_definition <- "standard_pando_coefficient"
      tab$analysis_mode <- "standard_pando"
      tab <- tab[, c(
        "group_id", condition_col, celltype_col,
        setdiff(colnames(tab), c("group_id", condition_col, celltype_col))
      ), drop = FALSE]
      tab
    }
    all_rows[[value]] <- add_design(extracted$all)
    active_rows[[value]] <- add_design(extracted$significant)
    status_rows[[value]] <- data.frame(
      group_id = rc_make_stratum_id(
        data.frame(
          condition_value = condition_levels[[1L]],
          celltype_value = value,
          stringsAsFactors = FALSE
        ),
        c("condition_value", "celltype_value")
      ),
      condition_value = condition_levels[[1L]],
      celltype_value = value,
      status = "ok",
      n_cells = length(cells),
      n_target_genes = length(target_genes),
      n_edges = nrow(extracted$all),
      n_active_edges = nrow(extracted$significant),
      grn_evidence_role = "standard_pando_full_fit",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    colnames(status_rows[[value]])[
      colnames(status_rows[[value]]) == "condition_value"
    ] <- condition_col
    colnames(status_rows[[value]])[
      colnames(status_rows[[value]]) == "celltype_value"
    ] <- celltype_col
    objects[[value]] <- grn
    selected_cells <- c(selected_cells, cells)
    if (isTRUE(save_pando_objects)) {
      saveRDS(grn, file.path(
        outdir, "pando_objects",
        paste0("standard_pando_", .rc_safe_path_component(value), ".rds")
      ))
    }
  }
  all_table <- .rc_bind_frames_fill(all_rows)
  active_table <- .rc_bind_frames_fill(active_rows)
  status <- do.call(rbind, status_rows)
  rownames(status) <- NULL
  .rc_step_monitor_event(
    progress_monitor, "standard_artifacts",
    "writing standard Pando Stage 1 artifacts", current = 10L,
    context = list(cell_types = length(cell_types))
  )
  .rc_write_tsv_gz(status, file.path(outdir, "pando_group_status.tsv.gz"))
  .rc_write_tsv_gz(all_table, file.path(outdir, "pando_tf_peak_gene_standard_all.tsv.gz"))
  .rc_write_tsv_gz(active_table, file.path(outdir, "pando_tf_peak_gene_standard_active.tsv.gz"))
  answer <- list(
    schema_version = "regcompass_standard_pando_fit_v1",
    analysis_mode = "standard_pando",
    condition_coefficients_calculated = FALSE,
    target_metabolic_genes = target_genes,
    condition_fit_status = status,
    standard_pando_objects = objects,
    condition_grn_fits = list(),
    pando_network_index = data.frame(),
    pando_fit_diagnostics = data.frame(),
    tf_peak_gene_universal = active_table,
    tf_peak_gene_condition_all = all_table,
    tf_peak_gene_condition = active_table,
    tf_peak_gene_condition_effect_all = data.frame(),
    tf_peak_gene_condition_effect = data.frame(),
    paired_cell_ids = selected_cells,
    paired_cell_metadata = data.frame(
      cell_id = selected_cells,
      condition = as.character(object@meta.data[selected_cells, condition_col]),
      cell_type = as.character(object@meta.data[selected_cells, celltype_col]),
      stringsAsFactors = FALSE
    ),
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF",
      grn_fit = "original Pando infer_grn per broad cell type",
      coefficient_contract = "standard_pando_coefficient_no_condition_effect",
      condition_coefficients_calculated = FALSE,
      standard_fallback_adjustments = fallback_adjustments,
      penalty_regulatory_evidence = paste(
        "standard Pando full-fit TF-by-ATAC cell projections aggregated by",
        "exact SuperCell membership"
      ),
      absolute_estimate_threshold = .rc_standard_pando_min_abs_fixed,
      padj_threshold = .rc_standard_pando_padj_fixed
    ),
    group_cols = c(condition_col, celltype_col)
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_standard_pando_projection <- function(
    grn_result, membership, metacell_meta, condition_col, celltype_col,
    rna_assay, atac_assay, target_genes) {
  units <- as.character(metacell_meta$unit_id %||% metacell_meta$metacell_id)
  projection <- matrix(
    NA_real_, length(target_genes), length(units),
    dimnames = list(tolower(target_genes), units)
  )
  reliability <- projection
  coverage <- list()
  for (celltype in names(grn_result$standard_pando_objects)) {
    grn <- grn_result$standard_pando_objects[[celltype]]
    cells <- colnames(grn@data)
    one_membership <- membership[
      membership$cell_id %in% cells, , drop = FALSE
    ]
    if (!nrow(one_membership)) next
    rna <- Matrix::t(Pando::LayerData(
      grn, assay = rna_assay, layer = "data"
    ))
    atac <- Matrix::t(Pando::LayerData(
      grn, assay = atac_assay, layer = "data"
    ))
    rownames(rna) <- colnames(grn@data)
    rownames(atac) <- colnames(grn@data)
    atac_keys <- .rc_pando_region_key(colnames(atac))
    if (anyDuplicated(atac_keys)) {
      stop("ATAC features are duplicated after Pando region normalization.",
           call. = FALSE)
    }
    colnames(atac) <- atac_keys
    edge <- grn_result$tf_peak_gene_condition[
      as.character(grn_result$tf_peak_gene_condition[[celltype_col]]) == celltype,
      , drop = FALSE
    ]
    if (!nrow(edge)) next
    edge$tf <- toupper(as.character(edge$tf))
    edge$target <- tolower(as.character(edge$target))
    edge$region_key <- .rc_pando_region_key(edge$region)
    rna_names <- toupper(colnames(rna))
    tf_index <- match(edge$tf, rna_names)
    peak_index <- match(edge$region_key, colnames(atac))
    usable <- !is.na(tf_index) & !is.na(peak_index) & is.finite(edge$estimate)
    edge <- edge[usable, , drop = FALSE]
    tf_index <- tf_index[usable]
    peak_index <- peak_index[usable]
    cell_score <- matrix(
      0, nrow = length(cells), ncol = length(target_genes),
      dimnames = list(cells, tolower(target_genes))
    )
    for (i in seq_len(nrow(edge))) {
      target <- edge$target[[i]]
      if (!target %in% colnames(cell_score)) next
      cell_score[, target] <- cell_score[, target] +
        as.numeric(rna[, tf_index[[i]]]) *
        as.numeric(atac[, peak_index[[i]]]) *
        as.numeric(edge$estimate[[i]])
    }
    mapped <- one_membership$cell_id %in% rownames(cell_score)
    one_membership <- one_membership[mapped, , drop = FALSE]
    groups <- split(one_membership$cell_id, one_membership$metacell_id)
    group_score <- do.call(rbind, lapply(groups, function(group_cells) {
      colMeans(cell_score[group_cells, , drop = FALSE])
    }))
    target <- intersect(colnames(group_score), rownames(projection))
    group <- intersect(rownames(group_score), colnames(projection))
    projection[target, group] <- t(group_score[group, target, drop = FALSE])
    rsq <- tapply(as.numeric(edge$rsq), edge$target, function(x) {
      x <- x[is.finite(x)]
      if (length(x)) max(x) else NA_real_
    })
    q <- sqrt(pmin(1, pmax(0, rsq[target])))
    reliability[target, group] <- matrix(
      q, nrow = length(target), ncol = length(group)
    )
    coverage[[celltype]] <- data.frame(
      target = target,
      cell_type = celltype,
      n_standard_pando_edges = vapply(target, function(gene) {
        sum(edge$target == gene)
      }, integer(1)),
      projection_origin = "standard_pando_full_fit",
      projection_used_for_penalty = TRUE,
      condition_coefficients_calculated = FALSE,
      stringsAsFactors = FALSE
    )
  }
  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    projection_origin = "standard_pando_full_fit",
    projection_used_for_penalty = TRUE,
    full_fit_projection_used_for_penalty = TRUE
  )
}
