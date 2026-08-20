# Merge independent condition-GRN cell-type results without merging GRNData feature spaces.

.rc_merge_condition_job_results <- function(results, full_object = NULL) {
  if (!is.list(results) || !length(results)) return(NULL)
  .rc_validate_pando_result_cell_partition(results)

  pando_objects <- lapply(results, function(value) value$pando_grn_data)
  valid_object <- vapply(
    pando_objects, function(value) inherits(value, "GRNData"), logical(1)
  )
  if (!all(valid_object)) {
    stop("Every condition-GRN cell-type result requires one Pando GRNData object.",
         call. = FALSE)
  }
  cell_types <- vapply(results, function(value) {
    fits <- value$condition_grn_fits
    labels <- unique(vapply(fits, function(fit) {
      as.character(fit$cell_type)[[1L]]
    }, character(1)))
    if (length(labels) != 1L || is.na(labels) || !nzchar(trimws(labels))) {
      stop("A condition-GRN result must contain exactly one cell type.",
           call. = FALSE)
    }
    labels
  }, character(1))
  if (anyDuplicated(cell_types)) {
    stop(
      "Condition-GRN cell-type results are duplicated: ",
      paste(unique(cell_types[duplicated(cell_types)]), collapse = ", "),
      call. = FALSE
    )
  }
  names(pando_objects) <- cell_types

  answer <- results[[1L]]
  if (length(results) > 1L) {
    frame_fields <- c(
      "condition_fit_status", "pando_network_index", "pando_fit_diagnostics",
      "pando_edge_inference", "tf_peak_gene_universal",
      "tf_peak_gene_condition_all", "tf_peak_gene_condition",
      "tf_peak_gene_condition_effect_all", "tf_peak_gene_condition_effect",
      "tf_peak_gene_condition_contrasts", "paired_cell_metadata"
    )
    for (field in frame_fields) {
      answer[[field]] <- .rc_bind_pando_field(results, field)
    }
    if (nrow(answer$paired_cell_metadata)) {
      answer$paired_cell_metadata <- answer$paired_cell_metadata[
        !duplicated(answer$paired_cell_metadata$cell_id), , drop = FALSE
      ]
    }
    answer$paired_cell_ids <- unique(
      as.character(answer$paired_cell_metadata$cell_id)
    )
    answer$target_metabolic_genes <- unique(unlist(
      lapply(results, `[[`, "target_metabolic_genes"), use.names = FALSE
    ))
    answer$condition_grn_fits <- unlist(
      lapply(results, `[[`, "condition_grn_fits"), recursive = FALSE
    )
    summaries <- lapply(results, `[[`, "pando_execution_summary")
    answer$pando_execution_summary <- list(
      fit_engine = paste(unique(unlist(lapply(
        summaries, `[[`, "fit_engine"
      ), use.names = FALSE)), collapse = ";"),
      targets_total = sum(vapply(summaries, function(x) {
        as.integer(x$targets_total %||% 0L)
      }, integer(1))),
      targets_failed = sum(vapply(summaries, function(x) {
        as.integer(x$targets_failed %||% 0L)
      }, integer(1)))
    )
  }

  answer$pando_grn_data_by_cell_type <- pando_objects
  answer$pando_grn_data <- if (length(pando_objects) == 1L) {
    pando_objects[[1L]]
  } else {
    NULL
  }
  answer$pando_object_scope <- list(
    cell_types = names(pando_objects),
    preserves_cell_type_peak_space = TRUE,
    combined_grndata = length(pando_objects) == 1L
  )

  fits <- answer$condition_grn_fits
  if (!is.list(fits) || !length(fits)) {
    stop("Merged condition-GRN output lacks fit contracts.", call. = FALSE)
  }
  fit_types <- vapply(fits, function(fit) as.character(fit$cell_type)[[1L]],
                      character(1))
  if (!setequal(fit_types, names(pando_objects))) {
    stop(
      "Merged condition fits and Pando objects cover different cell types.",
      call. = FALSE
    )
  }
  paired <- answer$paired_cell_metadata
  if (!is.data.frame(paired) || !"cell_id" %in% colnames(paired)) {
    stop("Merged condition-GRN output lacks paired-cell metadata.",
         call. = FALSE)
  }
  cell_id <- as.character(paired$cell_id)
  if (anyNA(cell_id) || any(!nzchar(trimws(cell_id))) || anyDuplicated(cell_id)) {
    stop("Merged condition-GRN output contains invalid paired-cell IDs.",
         call. = FALSE)
  }
  answer
}
