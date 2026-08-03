# Preserve exact per-cell-type Pando feature spaces after parallel Stage 1.

.rc_condition_job_object_map <- function(results) {
  if (!length(results)) return(list())
  values <- lapply(results, function(value) value$pando_grn_data)
  valid <- vapply(values, function(value) inherits(value, "GRNData"), logical(1))
  if (!all(valid)) {
    stop("Every condition-GRN job must return a Pando GRNData object.",
         call. = FALSE)
  }
  cell_types <- vapply(results, function(value) {
    fits <- value$condition_grn_fits
    labels <- unique(vapply(fits, function(fit) {
      as.character(fit$cell_type)[[1L]]
    }, character(1)))
    if (length(labels) != 1L || is.na(labels) || !nzchar(labels)) {
      stop("A condition-GRN job must contain exactly one cell type.",
           call. = FALSE)
    }
    labels
  }, character(1))
  if (anyDuplicated(cell_types)) {
    stop("Condition-GRN cell-type jobs are duplicated: ",
         paste(unique(cell_types[duplicated(cell_types)]), collapse = ", "),
         call. = FALSE)
  }
  names(values) <- cell_types
  values
}

.rc_merge_condition_job_results <- function(results, full_object = NULL) {
  if (!length(results)) return(NULL)
  answer <- results[[1L]]
  object_map <- .rc_condition_job_object_map(results)
  if (length(results) > 1L) {
    frame_fields <- c(
      "condition_fit_status", "pando_network_index", "pando_fit_diagnostics",
      "tf_peak_gene_universal", "tf_peak_gene_condition_all",
      "tf_peak_gene_condition", "tf_peak_gene_condition_effect_all",
      "tf_peak_gene_condition_effect", "paired_cell_metadata"
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
  answer$pando_grn_data_by_cell_type <- object_map
  answer$pando_grn_data <- if (length(object_map) == 1L) {
    object_map[[1L]]
  } else {
    NULL
  }
  answer$parallel_object_contract <- list(
    schema = "pando_grn_data_by_cell_type_v1",
    cell_types = names(object_map),
    preserves_local_peak_space = TRUE,
    merged_grndata_prohibited = length(object_map) > 1L
  )
  answer
}

.rc_merge_pando_results_base_parallel_contract <- .rc_merge_pando_results

.rc_merge_pando_results <- function(
    condition_result = NULL, standard_results = list(),
    condition_types = character(), standard_types = character(),
    condition_col, celltype_col, outdir) {
  answer <- .rc_merge_pando_results_base_parallel_contract(
    condition_result = condition_result,
    standard_results = standard_results,
    condition_types = condition_types,
    standard_types = standard_types,
    condition_col = condition_col,
    celltype_col = celltype_col,
    outdir = outdir
  )
  answer$pando_grn_data_by_cell_type <- if (is.null(condition_result)) {
    list()
  } else {
    condition_result$pando_grn_data_by_cell_type %||% list()
  }
  answer$parallel_object_contract <- if (is.null(condition_result)) {
    list(
      schema = "pando_grn_data_by_cell_type_v1",
      cell_types = character(),
      preserves_local_peak_space = TRUE,
      merged_grndata_prohibited = FALSE
    )
  } else {
    condition_result$parallel_object_contract
  }
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
