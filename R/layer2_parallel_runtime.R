# Layer 2 parallel context helpers.

.rc_layer2_parallel_context <- new.env(parent = emptyenv())
.rc_layer2_parallel_context$active <- FALSE
.rc_layer2_parallel_context$parallel <- TRUE
.rc_layer2_parallel_context$BPPARAM <- NULL
.rc_layer2_parallel_context$nested_serial <- FALSE

.rc_layer2_enter_parallel_context <- function(parallel, BPPARAM) {
  previous <- as.list(.rc_layer2_parallel_context)
  .rc_layer2_parallel_context$active <- TRUE
  .rc_layer2_parallel_context$parallel <- isTRUE(parallel)
  .rc_layer2_parallel_context$BPPARAM <- BPPARAM
  .rc_layer2_parallel_context$nested_serial <- FALSE
  previous
}

.rc_layer2_restore_parallel_context <- function(previous) {
  rm(
    list = ls(.rc_layer2_parallel_context, all.names = TRUE),
    envir = .rc_layer2_parallel_context
  )
  list2env(previous, envir = .rc_layer2_parallel_context)
  invisible(NULL)
}

.rc_atomic_save_rds <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(file), ".tmp_"),
    tmpdir = dirname(file)
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, file)) {
    if (!file.copy(temporary, file, overwrite = TRUE)) {
      stop("Cannot persist Layer 2 task result: ", file, call. = FALSE)
    }
    unlink(temporary, force = TRUE)
  }
  invisible(file)
}

.rc_layer2_task_bpparam <- function() {
  if (!isTRUE(.rc_layer2_parallel_context$active) ||
      !isTRUE(.rc_layer2_parallel_context$parallel) ||
      isTRUE(.rc_layer2_parallel_context$nested_serial)) {
    return(FALSE)
  }
  .rc_layer2_parallel_context$BPPARAM
}
