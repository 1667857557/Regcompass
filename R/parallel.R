# Parallel helpers used by the two-stage RegCompass workflow.

#' Detect a conservative RegCompass worker count
#' @export
rc_available_workers <- function(default = 1L) {
  vals <- c(
    getOption("RegCompassR.workers", NA),
    Sys.getenv("REGCOMPASS_WORKERS", unset = NA_character_),
    Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_),
    Sys.getenv("NSLOTS", unset = NA_character_)
  )
  vals <- suppressWarnings(as.integer(vals))
  vals <- vals[is.finite(vals) & vals >= 1L]
  if (length(vals)) return(max(1L, vals[[1L]]))
  if (requireNamespace("future", quietly = TRUE)) {
    value <- tryCatch(future::availableCores(), error = function(e) NA_integer_)
    value <- suppressWarnings(as.integer(value[[1L]]))
    if (is.finite(value) && value >= 1L) return(value)
  }
  value <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)[[1L]]))
  if (!is.finite(value) || value < 1L) max(1L, as.integer(default[[1L]])) else max(1L, value - 1L)
}

#' Build a RegCompass BiocParallel backend
#' @export
rc_default_bpparam <- function(workers = NULL, backend = c("auto", "serial", "snow", "multicore")) {
  backend <- match.arg(backend)
  if (is.null(workers)) workers <- rc_available_workers(default = 1L)
  workers <- suppressWarnings(as.integer(workers[[1L]]))
  if (is.na(workers) || workers < 2L || identical(backend, "serial")) return(NULL)
  if (!requireNamespace("BiocParallel", quietly = TRUE)) return(NULL)
  if (identical(backend, "auto")) {
    in_container <- nzchar(Sys.getenv("CONTAINER")) || file.exists("/.dockerenv")
    backend <- if (.Platform$OS.type == "windows" || in_container) "snow" else "multicore"
  }
  if (identical(backend, "snow")) {
    BiocParallel::SnowParam(workers = workers, type = "SOCK")
  } else {
    BiocParallel::MulticoreParam(workers = workers)
  }
}

.rc_stop_bpparam <- function(BPPARAM) {
  if (identical(BPPARAM, FALSE) || is.null(BPPARAM) ||
      !requireNamespace("BiocParallel", quietly = TRUE)) return(invisible(FALSE))
  started <- tryCatch(isTRUE(BiocParallel::bpstarted(BPPARAM)), error = function(e) FALSE)
  if (started) try(BiocParallel::bpstop(BPPARAM), silent = TRUE)
  invisible(started)
}

#' Apply a function with a stage-owned BiocParallel backend
#'
#' The backend is always stopped when the call returns. This creates a hard
#' process boundary between the upstream stratum stage and downstream LP stage.
#' @export
rc_parallel_lapply <- function(X, FUN, BPPARAM = NULL, ...) {
  if (identical(BPPARAM, FALSE) || length(X) <= 1L) return(lapply(X, FUN, ...))
  if (is.null(BPPARAM)) BPPARAM <- rc_default_bpparam()
  if (is.null(BPPARAM)) return(lapply(X, FUN, ...))
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel must be installed when `BPPARAM` is provided.", call. = FALSE)
  }
  if (!isTRUE(BiocParallel::bpstarted(BPPARAM))) BiocParallel::bpstart(BPPARAM)
  on.exit(.rc_stop_bpparam(BPPARAM), add = TRUE)
  BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)
}
