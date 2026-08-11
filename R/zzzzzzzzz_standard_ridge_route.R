# Optional standard-Pando ridge routing. Original GLM remains the default.

.rc_pando_infer_arg_catalog_standard_ridge_impl <- .rc_pando_infer_arg_catalog
.rc_pando_infer_arg_catalog <- function() {
  out <- .rc_pando_infer_arg_catalog_standard_ridge_impl()
  out$standard <- unique(c(
    out$standard,
    "ridge_control", "rank_action", "min_residual_df", "padj_threshold"
  ))
  out
}

.rc_route_pando_infer_args_standard_ridge_impl <- .rc_route_pando_infer_args
.rc_route_pando_infer_args <- function(
    args, condition_types = character(), standard_types = character()) {
  if (!is.list(args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  answer <- .rc_route_pando_infer_args_standard_ridge_impl(
    args,
    condition_types = condition_types,
    standard_types = standard_types
  )
  method <- as.character(answer$standard$method %||% "glm")
  if (length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("Standard Pando `method` must be one non-empty value.", call. = FALSE)
  }
  ridge_fields <- c("ridge_control", "rank_action", "min_residual_df",
                    "padj_threshold")
  shared_ridge_fields <- c("rank_action", "min_residual_df", "padj_threshold")

  if (identical(method, "ridge")) {
    answer$standard$ridge_control <- answer$standard$ridge_control %||% list()
    if (!is.list(answer$standard$ridge_control)) {
      stop("Standard Pando `ridge_control` must be a list.", call. = FALSE)
    }
    answer$standard$rank_action <- answer$standard$rank_action %||% "mark"
    answer$standard$min_residual_df <-
      answer$standard$min_residual_df %||% 1L
    answer$standard$padj_threshold <-
      answer$standard$padj_threshold %||% 0.05
  } else {
    if ("ridge_control" %in% names(args) && length(args$ridge_control)) {
      stop(
        "`pando_infer_args$ridge_control` requires standard Pando ",
        "`method = \"ridge\"`.", call. = FALSE
      )
    }
    answer$standard[ridge_fields] <- NULL
  }

  if (is.data.frame(answer$diagnostics) && nrow(answer$diagnostics)) {
    false_condition_disable <-
      answer$diagnostics$route == "condition_grn" &
      answer$diagnostics$argument %in% shared_ridge_fields
    valid_standard_ridge <- identical(method, "ridge") &
      answer$diagnostics$route == "standard_pando" &
      answer$diagnostics$argument %in% ridge_fields
    answer$diagnostics <- answer$diagnostics[
      !(false_condition_disable | valid_standard_ridge),
      , drop = FALSE
    ]
    rownames(answer$diagnostics) <- NULL
  }
  answer
}

.rc_standard_pando_infer_args_ridge_impl <- .rc_standard_pando_infer_args
.rc_standard_pando_infer_args <- function(args) {
  if (!is.list(args)) stop("`pando_infer_args` must be a list.", call. = FALSE)
  method <- as.character(args$method %||% "glm")
  is_ridge <- length(method) == 1L && !is.na(method) && identical(method, "ridge")
  target_param <- args$BPPARAM %||% NULL
  if (is_ridge) args$BPPARAM <- NULL
  answer <- .rc_standard_pando_infer_args_ridge_impl(args)
  if (is_ridge && !is.null(target_param) && !identical(target_param, FALSE)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE) ||
        !methods::is(target_param, "BiocParallelParam")) {
      stop("Standard ridge `BPPARAM` must be a BiocParallelParam.",
           call. = FALSE)
    }
    answer$BPPARAM <- target_param
  }
  answer
}

.rc_run_standard_pando_celltype_job_glm_impl <-
  .rc_run_standard_pando_celltype_job

.rc_run_standard_pando_celltype_job <- function(
    job, base, extra_args, standard_infer_args,
    outer_parallel, progress_monitor) {
  method <- as.character(standard_infer_args$method %||% "glm")
  if (!identical(method, "ridge")) {
    return(.rc_run_standard_pando_celltype_job_glm_impl(
      job = job,
      base = base,
      extra_args = extra_args,
      standard_infer_args = standard_infer_args,
      outer_parallel = outer_parallel,
      progress_monitor = progress_monitor
    ))
  }
  if (!is.list(job) || !inherits(job$object, "Seurat")) {
    stop("Invalid standard-Pando cell-type job.", call. = FALSE)
  }

  thread_state <- .rc_set_internal_single_thread()
  target_param <- NULL
  on.exit({
    if (!is.null(target_param) &&
        requireNamespace("BiocParallel", quietly = TRUE) &&
        isTRUE(tryCatch(BiocParallel::bpisup(target_param),
                        error = function(e) FALSE))) {
      try(BiocParallel::bpstop(target_param), silent = TRUE)
    }
    .rc_restore_internal_threads(thread_state)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  grn_budget <- suppressWarnings(as.integer(
    Sys.getenv("REGCOMPASS_TASK_WORKERS", unset = "1")
  ))
  if (!is.finite(grn_budget) || grn_budget < 1L) grn_budget <- 1L
  target_workers <- if (isTRUE(outer_parallel)) {
    max(0L, grn_budget - 1L)
  } else {
    grn_budget
  }
  target_parallel <- target_workers >= 2L
  if (target_parallel) {
    target_param <- .rc_condition_nested_target_bpparam(target_workers)
  }

  job_extra <- extra_args
  motif_args <- job_extra$pando_motif_args %||% list()
  if (outer_parallel && is.list(motif_args) &&
      !is.null(motif_args$cache_dir)) {
    motif_args$cache_dir <- file.path(
      motif_args$cache_dir, .rc_safe_path_component(job$cell_type)
    )
    job_extra$pando_motif_args <- motif_args
  }
  args <- c(base[setdiff(names(base), names(job_extra))], job_extra)
  args$object <- job$object
  args$cell_type <- job$cell_type
  args$progress_monitor <- if (outer_parallel) NULL else progress_monitor
  args$outdir <- file.path(
    base$outdir, "standard", .rc_safe_path_component(job$cell_type)
  )

  fixed_infer_args <- .rc_standard_pando_single_process_args(
    standard_infer_args
  )
  original_contract <- attr(
    fixed_infer_args, "regcompass_pando_parallel_contract", exact = TRUE
  )
  attr(fixed_infer_args, "regcompass_pando_parallel_contract") <- NULL
  if (target_parallel) fixed_infer_args$BPPARAM <- target_param
  args$pando_infer_args <- fixed_infer_args
  args$parallel <- target_parallel

  value <- do.call(.rc_fit_standard_pando_by_cell_type, args)
  parallel_contract <- utils::modifyList(
    original_contract %||% list(),
    list(
      scope = "standard_pando_ridge_target",
      method = "ridge",
      grn_worker_budget = as.integer(grn_budget),
      outer_celltype_worker = isTRUE(outer_parallel),
      target_workers = as.integer(if (target_parallel) target_workers else 1L),
      target_parallel = target_parallel,
      target_backend = if (target_parallel) "snow" else "serial",
      target_pool_released_after_cell_type = TRUE,
      worker_budget_bounded = TRUE
    )
  )
  if (is.list(value$normalization_policy)) {
    value$normalization_policy$parallel_contract <- parallel_contract
  }
  list(cell_type = job$cell_type, route = "standard_pando", result = value)
}
