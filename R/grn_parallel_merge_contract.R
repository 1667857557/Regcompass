# Strict contract for merging independent condition-GRN worker results.

.rc_merge_condition_job_results_contract_base <-
  .rc_merge_condition_job_results

.rc_validate_condition_job_merge_output <- function(answer) {
  if (is.null(answer)) return(invisible(TRUE))
  fits <- answer$condition_grn_fits
  object_map <- answer$pando_grn_data_by_cell_type
  if (!is.list(fits) || !length(fits) ||
      !is.list(object_map) || !length(object_map)) {
    stop("Parallel condition-GRN output lacks fits or per-cell-type objects.",
         call. = FALSE)
  }
  fit_types <- vapply(fits, function(fit) {
    value <- as.character(fit$cell_type)
    if (length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
      stop("A parallel condition fit lacks one cell-type label.",
           call. = FALSE)
    }
    value
  }, character(1))
  map_names <- names(object_map)
  if (anyDuplicated(fit_types) || is.null(map_names) ||
      anyNA(map_names) || any(!nzchar(map_names)) ||
      anyDuplicated(map_names) || !setequal(fit_types, map_names)) {
    stop(
      "Parallel condition fits and per-cell-type Pando objects do not form ",
      "the same unique cell-type set.", call. = FALSE
    )
  }
  invalid_object <- vapply(object_map, function(value) {
    !inherits(value, "GRNData")
  }, logical(1))
  if (any(invalid_object)) {
    stop("A parallel condition cell type lacks its Pando GRNData object.",
         call. = FALSE)
  }
  paired <- answer$paired_cell_metadata
  if (!is.data.frame(paired) || !"cell_id" %in% colnames(paired)) {
    stop("Parallel condition output lacks paired-cell metadata.",
         call. = FALSE)
  }
  cell_id <- as.character(paired$cell_id)
  if (anyNA(cell_id) || any(!nzchar(trimws(cell_id))) ||
      anyDuplicated(cell_id)) {
    stop("Parallel condition output contains invalid or duplicated cell IDs.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_merge_condition_job_results <- function(results, full_object = NULL) {
  .rc_validate_pando_result_cell_partition(results)
  answer <- .rc_merge_condition_job_results_contract_base(
    results = results, full_object = full_object
  )
  .rc_validate_condition_job_merge_output(answer)
  answer
}
