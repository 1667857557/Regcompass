# Runtime and schema hardening for the shared-backbone GRN implementation.

.rc_mt_validate_args_core <- .rc_mt_validate_args
.rc_mt_validate_args <- function(args) {
  answer <- .rc_mt_validate_args_core(args)
  if (answer$min_cv_rsq > 1) {
    stop("`multitask_args$min_cv_rsq` must be in [0, 1].", call. = FALSE)
  }
  answer
}

.rc_mt_empty_target_fit <- function(target_design, condition, reason, message) {
  edge <- target_design$edge_metadata
  if (!is.data.frame(edge)) edge <- data.frame()
  conditions <- target_design$condition_levels %||%
    sort(unique(as.character(condition)))
  edge_ids <- as.character(edge$edge_id %||% character())
  zero <- matrix(
    0,
    nrow = length(edge_ids),
    ncol = length(conditions),
    dimnames = list(edge_ids, conditions)
  )
  edge$fit_status <- rep(reason, nrow(edge))
  edge$fit_message <- rep(as.character(message), nrow(edge))
  list(
    effects = list(
      global = stats::setNames(rep(0, length(edge_ids)), edge_ids),
      deviation = zero,
      effective = zero,
      internal_global = stats::setNames(rep(0, length(edge_ids)), edge_ids),
      internal_deviation = zero
    ),
    selection_frequency = zero,
    sign_stability = zero,
    lambda = NA_real_,
    cv_rsq = NA_real_,
    n_stability_successful = 0L,
    n_stability_requested = 0L,
    edge_metadata = edge,
    fit_status = reason,
    fit_message = as.character(message)
  )
}

.rc_mt_fit_target_core <- .rc_mt_fit_target
.rc_mt_fit_target <- function(y, target_design, condition, args, target) {
  y <- as.numeric(y)
  valid_y <- y[is.finite(y)]
  if (!isTRUE(target_design$estimable)) {
    return(.rc_mt_empty_target_fit(
      target_design, condition,
      reason = "skipped_no_estimable_interaction",
      message = paste0("No within-condition interaction variance for ", target)
    ))
  }
  if (length(valid_y) < 3L || !is.finite(stats::var(valid_y)) ||
      stats::var(valid_y) <= args$coefficient_zero_tolerance) {
    return(.rc_mt_empty_target_fit(
      target_design, condition,
      reason = "skipped_zero_target_variance",
      message = paste0("Target expression has no estimable variance for ", target)
    ))
  }
  answer <- tryCatch(
    .rc_mt_fit_target_core(y, target_design, condition, args, target),
    error = function(error) error
  )
  if (inherits(answer, "error")) {
    return(.rc_mt_empty_target_fit(
      target_design, condition,
      reason = "skipped_model_failure",
      message = conditionMessage(answer)
    ))
  }
  answer$edge_metadata$fit_status <- rep("ok", nrow(answer$edge_metadata))
  answer$edge_metadata$fit_message <- rep(NA_character_, nrow(answer$edge_metadata))
  answer$fit_status <- "ok"
  answer$fit_message <- NA_character_
  answer
}

.rc_mt_effective_direction <- function(x, tolerance) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    !is.finite(x) | abs(x) <= tolerance,
    "zero",
    ifelse(x > 0, "positive", "negative")
  )
}

.rc_mt_mark_sign_flips <- function(edges, celltype_col, tolerance) {
  if (!is.data.frame(edges) || !nrow(edges)) return(edges)
  required <- c(celltype_col, "edge_id", "effective_estimate", "active_edge")
  if (!all(required %in% colnames(edges))) return(edges)
  key <- paste(
    as.character(edges[[celltype_col]]),
    as.character(edges$edge_id),
    sep = "\001"
  )
  rows <- split(seq_len(nrow(edges)), key)
  flag <- rep(FALSE, nrow(edges))
  for (index in rows) {
    selected <- index[edges$active_edge[index] %in% TRUE]
    value <- suppressWarnings(as.numeric(edges$effective_estimate[selected]))
    value <- value[is.finite(value) & abs(value) > tolerance]
    flip <- length(value) >= 2L && any(value > 0) && any(value < 0)
    flag[index] <- flip
  }
  edges$sign_flip_flag <- flag
  edges$effective_direction <- .rc_mt_effective_direction(
    edges$effective_estimate, tolerance
  )
  edges
}

.rc_mt_condition_target_table <- function(
    significant, condition_col, celltype_col) {
  columns <- c("group_id", condition_col, celltype_col, "target")
  if (!is.data.frame(significant) || !nrow(significant) ||
      !all(columns %in% colnames(significant))) {
    answer <- data.frame(
      group_id = character(),
      condition = character(),
      cell_type = character(),
      target = character(),
      n_active_edges = integer(),
      n_positive_edges = integer(),
      n_negative_edges = integer(),
      stringsAsFactors = FALSE
    )
    names(answer)[2:3] <- c(condition_col, celltype_col)
    return(answer)
  }
  key <- paste(
    significant$group_id,
    as.character(significant$target),
    sep = "\001"
  )
  rows <- split(seq_len(nrow(significant)), key)
  answer <- do.call(rbind, lapply(rows, function(index) {
    one <- significant[index, , drop = FALSE]
    estimate <- suppressWarnings(as.numeric(one$effective_estimate))
    values <- one[1L, c("group_id", condition_col, celltype_col), drop = FALSE]
    data.frame(
      values,
      target = as.character(one$target[[1L]]),
      n_active_edges = nrow(one),
      n_positive_edges = sum(is.finite(estimate) & estimate > 0),
      n_negative_edges = sum(is.finite(estimate) & estimate < 0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  rownames(answer) <- NULL
  answer
}

.rc_run_celltype_multitask_grns_core <- .rc_run_celltype_multitask_grns
.rc_run_celltype_multitask_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    sample_col = NULL,
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_design_args = list(
      screen_method = "structural",
      min_tf_detection = 0.01,
      min_peak_detection = 0.01,
      min_target_detection = 0.01
    ),
    multitask_args = list(),
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    species = c("auto", "human", "mouse")) {
  meta <- object@meta.data
  conditions <- sort(unique(trimws(as.character(meta[[condition_col]]))))
  if (length(conditions) < 2L) {
    stop(
      "`multitask_shared_backbone` requires at least two conditions.",
      call. = FALSE
    )
  }
  cell_types <- sort(unique(trimws(as.character(meta[[celltype_col]]))))
  for (cell_type in cell_types) {
    observed <- sort(unique(trimws(as.character(
      meta[[condition_col]][as.character(meta[[celltype_col]]) == cell_type]
    ))))
    missing <- setdiff(conditions, observed)
    if (length(missing)) {
      stop(
        "Cell type ", cell_type,
        " is absent from conditions required by the shared task design: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (!is.null(sample_col)) {
    sample <- trimws(as.character(meta[[sample_col]]))
    if (anyNA(sample) || any(!nzchar(sample))) {
      stop("Sample metadata must be complete in shared-backbone mode.",
           call. = FALSE)
    }
  }

  answer <- .rc_run_celltype_multitask_grns_core(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    condition_col = condition_col,
    celltype_col = celltype_col,
    sample_col = sample_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    min_cells = min_cells,
    pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args,
    pando_design_args = pando_design_args,
    multitask_args = multitask_args,
    save_pando_objects = save_pando_objects,
    BPPARAM = BPPARAM,
    species = species
  )
  args <- answer$multitask_model$args
  design_map <- unique(answer$tf_peak_gene_candidates[, c(
    celltype_col, "design_id"
  ), drop = FALSE])
  attach_design <- function(data) {
    if (!is.data.frame(data) || !nrow(data) ||
        !celltype_col %in% colnames(data)) return(data)
    data$design_id <- design_map$design_id[
      match(as.character(data[[celltype_col]]),
            as.character(design_map[[celltype_col]]))
    ]
    data
  }
  answer$tf_peak_gene_global <- attach_design(answer$tf_peak_gene_global)
  answer$tf_peak_gene_condition_all <- attach_design(
    answer$tf_peak_gene_condition_all
  )
  answer$tf_peak_gene_all <- answer$tf_peak_gene_condition_all
  answer$tf_peak_gene_significant <- attach_design(
    answer$tf_peak_gene_significant
  )
  answer$tf_peak_gene_condition_all <- .rc_mt_mark_sign_flips(
    answer$tf_peak_gene_condition_all,
    celltype_col = celltype_col,
    tolerance = args$coefficient_zero_tolerance
  )
  answer$tf_peak_gene_all <- answer$tf_peak_gene_condition_all
  if (nrow(answer$tf_peak_gene_significant)) {
    flip_lookup <- unique(answer$tf_peak_gene_condition_all[, c(
      celltype_col, "edge_id", "sign_flip_flag"
    ), drop = FALSE])
    flip_key <- paste(
      as.character(flip_lookup[[celltype_col]]), flip_lookup$edge_id,
      sep = "\001"
    )
    significant_key <- paste(
      as.character(answer$tf_peak_gene_significant[[celltype_col]]),
      answer$tf_peak_gene_significant$edge_id,
      sep = "\001"
    )
    answer$tf_peak_gene_significant$sign_flip_flag <-
      flip_lookup$sign_flip_flag[match(significant_key, flip_key)]
    answer$tf_peak_gene_significant$effective_direction <-
      .rc_mt_effective_direction(
        answer$tf_peak_gene_significant$effective_estimate,
        args$coefficient_zero_tolerance
      )
  }
  if (!nrow(answer$tf_peak_gene_significant)) {
    stop(
      "No stability-selected condition TF-peak-gene edges passed the ",
      "multi-task thresholds.",
      call. = FALSE
    )
  }
  answer$condition_target_genes <- .rc_mt_condition_target_table(
    answer$tf_peak_gene_significant, condition_col, celltype_col
  )

  status <- answer$celltype_fit_status
  if (is.data.frame(status) && nrow(status)) {
    all_edges <- answer$tf_peak_gene_condition_all
    fit_rows <- split(
      seq_len(nrow(all_edges)),
      paste(as.character(all_edges[[celltype_col]]), all_edges$target,
            sep = "\001")
    )
    fit_lookup <- lapply(fit_rows, function(index) {
      one <- all_edges[index, , drop = FALSE]
      list(
        status = as.character(one$fit_status[[1L]] %||% "ok"),
        message = as.character(one$fit_message[[1L]] %||% NA_character_)
      )
    })
    for (i in seq_len(nrow(status))) {
      key <- paste(
        as.character(status[[celltype_col]][[i]]),
        as.character(status$target[[i]]),
        sep = "\001"
      )
      info <- fit_lookup[[key]]
      if (!is.null(info)) {
        status$status[[i]] <- info$status
        status$error_message[[i]] <- info$message
      }
    }
    answer$celltype_fit_status <- status
  }

  .rc_mm_write_tsv_gz(
    answer$celltype_fit_status,
    file.path(outdir, "multitask_target_status.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    answer$tf_peak_gene_global,
    file.path(outdir, "pando_tf_peak_gene_global.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    answer$tf_peak_gene_condition_all,
    file.path(outdir, "pando_tf_peak_gene_all.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    answer$tf_peak_gene_significant,
    file.path(outdir, "pando_tf_peak_gene_significant.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    answer$condition_target_genes,
    file.path(outdir, "condition_target_genes.tsv.gz")
  )
  saveRDS(answer, file.path(outdir, "multitask_grn_result.rds"))
  answer
}
