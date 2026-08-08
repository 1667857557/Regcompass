# Validation helpers for condition/standard Pando routing.

.rc_validate_pando_route_partition <- function(
    condition_types = character(), standard_types = character()) {
  normalize <- function(value, label) {
    value <- as.character(value)
    if (anyNA(value) || any(!nzchar(trimws(value))) || anyDuplicated(value)) {
      stop(label, " cell-type routes must be complete and unique.",
           call. = FALSE)
    }
    value
  }
  condition_types <- normalize(condition_types, "Condition-Pando")
  standard_types <- normalize(standard_types, "Standard-Pando")
  overlap <- intersect(condition_types, standard_types)
  if (length(overlap)) {
    stop(
      "Cell types cannot be assigned to both condition and standard Pando: ",
      paste(overlap, collapse = ", "), ".", call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_validate_pando_result_cell_partition <- function(results) {
  if (!is.list(results) || !length(results)) return(invisible(TRUE))
  metadata <- lapply(results, function(value) value$paired_cell_metadata)
  metadata <- metadata[vapply(metadata, is.data.frame, logical(1))]
  metadata <- metadata[vapply(metadata, nrow, integer(1)) > 0L]
  if (!length(metadata)) return(invisible(TRUE))
  missing_id <- vapply(metadata, function(value) {
    !"cell_id" %in% colnames(value)
  }, logical(1))
  if (any(missing_id)) {
    stop("A Pando route result lacks paired-cell IDs.", call. = FALSE)
  }
  cell_id <- as.character(unlist(
    lapply(metadata, `[[`, "cell_id"), use.names = FALSE
  ))
  if (anyNA(cell_id) || any(!nzchar(trimws(cell_id)))) {
    stop("Pando route cell IDs must be complete non-empty strings.",
         call. = FALSE)
  }
  duplicated_id <- unique(cell_id[duplicated(cell_id)])
  if (length(duplicated_id)) {
    stop(
      "A paired cell was emitted by more than one Pando route; first duplicate: ",
      duplicated_id[[1L]], ".", call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_overlay_projection <- function(target, incoming, label) {
  if (!identical(dimnames(target), dimnames(incoming))) {
    stop(label, " projection layout is incompatible.", call. = FALSE)
  }
  overlap <- is.finite(target) & is.finite(incoming)
  if (any(overlap)) {
    index <- which(overlap, arr.ind = TRUE)[1L, , drop = TRUE]
    stop(
      label, " projection overlaps an existing route at gene '",
      rownames(target)[index[[1L]]], "' and unit '",
      colnames(target)[index[[2L]]], "'.", call. = FALSE
    )
  }
  supplied <- is.finite(incoming)
  target[supplied] <- incoming[supplied]
  target
}
