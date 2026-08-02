# Late-loaded workflow hardening overrides. This file is collated last so the
# public APIs retain their names while replacing the earlier implementations.

.rc_stage1_min_cells_fixed <- 300L
.rc_standard_pando_padj_fixed <- 0.05
.rc_standard_pando_min_abs_fixed <- 0.01
.rc_lp_target_block_default <- 32L

.rc_original_step_grn_hardening <- rc_regcompass_step_grn
.rc_original_standard_pando_hardening <- .rc_fit_standard_pando_by_cell_type
.rc_original_condition_pando_hardening <- .rc_fit_condition_grns_by_cell_type
.rc_original_step_layer2_hardening <- rc_regcompass_step_layer2
.rc_original_run_microcompass_hardening <- rc_run_microcompass
.rc_original_run_microcompass_engine_hardening <- .rc_run_microcompass_engine
.rc_original_union_cache_hardening <- .rc_build_medium_specific_union_gem_cache
.rc_original_step_results_hardening <- rc_regcompass_step_results
.rc_original_extract_pando_hardening <- rc_extract_pando_tf_peak_gene

.rc_layer2_runtime_context <- new.env(parent = emptyenv())
.rc_layer2_runtime_context$monitor <- NULL
.rc_layer2_runtime_context$parallel <- FALSE
.rc_layer2_runtime_context$BPPARAM <- FALSE
.rc_layer2_runtime_context$directional_call <- 0L

.rc_prefilter_stage1_celltypes <- function(
    object, celltype_col, cell_type = NULL,
    min_cells = .rc_stage1_min_cells_fixed) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  .rc_validate_celltype_metadata(object@meta.data, celltype_col)
  observed <- trimws(as.character(object@meta.data[[celltype_col]]))
  available <- unique(observed)
  requested <- if (is.null(cell_type)) available else unique(trimws(as.character(cell_type)))
  missing <- setdiff(requested, available)
  if (length(missing)) {
    stop(
      "Requested cell types were not found: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  counts <- table(factor(observed[observed %in% requested], levels = requested))
  keep_types <- names(counts)[as.integer(counts) >= as.integer(min_cells)]
  drop_types <- names(counts)[as.integer(counts) < as.integer(min_cells)]
  if (!length(keep_types)) {
    detail <- paste0(names(counts), "=", as.integer(counts), collapse = ", ")
    stop(
      "No requested cell type reaches the fixed Stage 1 threshold of ",
      min_cells, " cells. Observed: ", detail,
      call. = FALSE
    )
  }
  if (length(drop_types)) {
    message(
      "Stage 1 excluded cell types below 300 cells before normalization: ",
      paste0(drop_types, "=", as.integer(counts[drop_types]), collapse = ", ")
    )
  }
  keep_cells <- rownames(object@meta.data)[observed %in% keep_types]
  filtered <- subset(object, cells = keep_cells)
  diagnostics <- data.frame(
    cell_type = names(counts),
    n_cells = as.integer(counts),
    retained = names(counts) %in% keep_types,
    threshold = as.integer(min_cells),
    stringsAsFactors = FALSE
  )
  filtered@misc$regcompass_stage1_group_filter <- diagnostics
  list(
    object = filtered,
    retained_cell_types = keep_types,
    diagnostics = diagnostics
  )
}

.rc_validate_stage1_fragment_policy <- function(fragment_files) {
  if (is.null(fragment_files)) return(FALSE)
  if (!is.logical(fragment_files) || length(fragment_files) != 1L ||
      is.na(fragment_files)) {
    stop(
      "Stage 1 `fragment_files` must be TRUE or FALSE; FALSE clears stale Signac fragment references.",
      call. = FALSE
    )
  }
  isTRUE(fragment_files)
}

.rc_pando_supports_motif_cache <- function() {
  if (!requireNamespace("Pando", quietly = TRUE)) return(FALSE)
  method <- tryCatch(
    getS3method("find_motifs", "GRNData", envir = asNamespace("Pando")),
    error = function(error) NULL
  )
  is.function(method) && all(c("cache_dir", "reuse_cache") %in% names(formals(method)))
}

#' Infer regulatory evidence after fixed Stage 1 prefiltering
#'
#' Cell types with fewer than 300 paired cells are removed before RNA/ATAC
#' normalization. Globally zero ATAC peaks are removed before TF-IDF or Pando.
#' By default stale Signac fragment references are cleared because Stage 1 uses
#' the in-memory peak matrix and genome sequence, not fragment files.
#'
#' @param fragment_files Preserve existing Signac fragment references when TRUE.
#' The default FALSE clears them before Stage 1.
#' @export
rc_regcompass_step_grn <- function(
    object, gem, outdir, genome,
    pfm = NULL,
    species = c("auto", "human", "mouse"),
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    pando_args = list(),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  if (!is.list(pando_args)) stop("`pando_args` must be a list.", call. = FALSE)
  preserve_fragments <- .rc_validate_stage1_fragment_policy(fragment_files)
  prefiltered <- .rc_prefilter_stage1_celltypes(
    object = object,
    celltype_col = celltype_col,
    cell_type = cell_type,
    min_cells = .rc_stage1_min_cells_fixed
  )
  object <- prefiltered$object
  if (!preserve_fragments) {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  }
  zero_filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Stage 1 prefilter"
  )
  object <- zero_filtered$object
  object@misc$regcompass_stage1_zero_peak_filter <- zero_filtered$diagnostics
  object@misc$regcompass_stage1_fragment_policy <- list(
    fragment_files = preserve_fragments,
    policy = if (preserve_fragments) "preserve" else "clear_before_stage1"
  )

  supplied_min <- pando_args$min_cells %||% .rc_stage1_min_cells_fixed
  supplied_min <- suppressWarnings(as.integer(supplied_min[[1L]]))
  if (!is.finite(supplied_min) || supplied_min != .rc_stage1_min_cells_fixed) {
    message("Stage 1 `min_cells` is fixed at 300; overriding the supplied value.")
  }
  pando_args$min_cells <- .rc_stage1_min_cells_fixed
  motif_args <- pando_args$pando_motif_args %||% list()
  if (!is.list(motif_args)) {
    stop("`pando_motif_args` must be a list.", call. = FALSE)
  }
  if (.rc_pando_supports_motif_cache()) {
    motif_args$cache_dir <- motif_args$cache_dir %||%
      file.path(outdir, "motif_cache")
    motif_args$reuse_cache <- motif_args$reuse_cache %||% TRUE
  }
  pando_args$pando_motif_args <- motif_args

  duplicate_file <- file.path(outdir, "single_cell_grn.rds")
  if (file.exists(duplicate_file)) unlink(duplicate_file, force = TRUE)
  on.exit({
    if (file.exists(duplicate_file)) unlink(duplicate_file, force = TRUE)
  }, add = TRUE)

  .rc_original_step_grn_hardening(
    object = object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = if (is.null(cell_type)) NULL else prefiltered$retained_cell_types,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    pando_args = pando_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
}

.rc_skip_duplicate_grn_save <- function(object, file, ...) {
  if (identical(basename(as.character(file)), "single_cell_grn.rds")) {
    return(invisible(NULL))
  }
  base::saveRDS(object, file, ...)
}

.rc_standard_pando_extract_strict <- function(...) {
  args <- list(...)
  requested_abs <- suppressWarnings(as.numeric(
    args$min_abs_estimate %||% .rc_standard_pando_min_abs_fixed
  ))
  if (!is.finite(requested_abs)) requested_abs <- .rc_standard_pando_min_abs_fixed
  abs_threshold <- max(.rc_standard_pando_min_abs_fixed, requested_abs)
  min_rsq <- suppressWarnings(as.numeric(args$min_model_rsq %||% 0))
  if (!is.finite(min_rsq)) min_rsq <- 0
  args$padj_threshold <- .rc_standard_pando_padj_fixed
  args$min_abs_estimate <- 0
  args$require_padj <- TRUE
  answer <- do.call(.rc_original_extract_pando_hardening, args)
  tab <- answer$all
  required <- c("estimate", "padj", "rsq")
  if (!all(required %in% colnames(tab))) {
    stop(
      "Standard Pando requires estimate, padj and rsq for strict edge filtering.",
      call. = FALSE
    )
  }
  estimate <- suppressWarnings(as.numeric(tab$estimate))
  padj <- suppressWarnings(as.numeric(tab$padj))
  rsq <- suppressWarnings(as.numeric(tab$rsq))
  keep <- is.finite(estimate) & abs(estimate) > abs_threshold &
    is.finite(padj) & padj < .rc_standard_pando_padj_fixed &
    is.finite(rsq) & rsq >= min_rsq
  answer$significant <- tab[keep, , drop = FALSE]
  attr(answer$significant, "edge_filter") <- list(
    padj = "< 0.05",
    min_abs_estimate = paste0("> ", format(abs_threshold, scientific = FALSE)),
    min_model_rsq = min_rsq
  )
  answer
}

.rc_fit_standard_pando_by_cell_type <- function(...) {
  args <- list(...)
  requested_abs <- suppressWarnings(as.numeric(
    args$min_abs_estimate %||% .rc_standard_pando_min_abs_fixed
  ))
  if (!is.finite(requested_abs)) requested_abs <- .rc_standard_pando_min_abs_fixed
  args$min_abs_estimate <- max(.rc_standard_pando_min_abs_fixed, requested_abs)
  implementation <- .rc_original_standard_pando_hardening
  evaluation_environment <- new.env(parent = environment(implementation))
  evaluation_environment$saveRDS <- .rc_skip_duplicate_grn_save
  evaluation_environment$rc_extract_pando_tf_peak_gene <-
    .rc_standard_pando_extract_strict
  environment(implementation) <- evaluation_environment
  do.call(implementation, args)
}

.rc_fit_condition_grns_by_cell_type <- function(...) {
  implementation <- .rc_original_condition_pando_hardening
  evaluation_environment <- new.env(parent = environment(implementation))
  evaluation_environment$saveRDS <- .rc_skip_duplicate_grn_save
  environment(implementation) <- evaluation_environment
  implementation(...)
}

# Rebuild the TF-IDF matrix once from triplets instead of repeatedly assigning
# non-contiguous row blocks into a large dgCMatrix.
.rc_apply_celltype_shared_tfidf <- function(
    object, celltype_col, atac_assay = "ATAC",
    method = 1, scale.factor = 1e4) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  .rc_validate_celltype_metadata(object@meta.data, celltype_col)
  object <- .rc_prepare_seurat_assays(
    object,
    assays = atac_assay,
    required_layers = "counts",
    optional_layers = "data"
  )
  filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Cell-type-shared TF-IDF"
  )
  object <- filtered$object
  counts <- .rc_get_assay_counts(object, atac_assay)
  units <- colnames(counts)
  meta <- object@meta.data[
    match(units, rownames(object@meta.data)), , drop = FALSE
  ]
  cell_type <- trimws(as.character(meta[[celltype_col]]))
  groups <- split(units, cell_type)
  local_zero <- vapply(groups, function(group_units) {
    sum(Matrix::rowSums(counts[, group_units, drop = FALSE]) <= 0)
  }, integer(1))
  row_index <- vector("list", length(groups))
  col_index <- vector("list", length(groups))
  values <- vector("list", length(groups))
  names(row_index) <- names(col_index) <- names(values) <- names(groups)
  for (group_name in names(groups)) {
    group_units <- groups[[group_name]]
    group_counts <- counts[, group_units, drop = FALSE]
    keep <- Matrix::rowSums(group_counts) > 0
    if (!any(keep)) next
    normalized <- Signac::RunTFIDF(
      group_counts[keep, , drop = FALSE],
      method = method,
      scale.factor = scale.factor,
      verbose = FALSE
    )
    triplet <- methods::as(.rc_as_sparse(normalized), "TsparseMatrix")
    if (!length(triplet@x)) next
    row_index[[group_name]] <- which(keep)[triplet@i + 1L]
    group_columns <- match(group_units, units)
    col_index[[group_name]] <- group_columns[triplet@j + 1L]
    values[[group_name]] <- as.numeric(triplet@x)
  }
  i <- unlist(row_index, use.names = FALSE)
  j <- unlist(col_index, use.names = FALSE)
  x <- unlist(values, use.names = FALSE)
  if (length(x)) {
    tfidf <- Matrix::sparseMatrix(
      i = i, j = j, x = x,
      dims = dim(counts),
      dimnames = dimnames(counts),
      index1 = TRUE,
      giveCsparse = TRUE
    )
  } else {
    tfidf <- Matrix::Matrix(
      0,
      nrow = nrow(counts),
      ncol = ncol(counts),
      sparse = TRUE,
      dimnames = dimnames(counts)
    )
  }
  tfidf <- .rc_as_sparse(tfidf)
  if (!.rc_same_matrix_layout(tfidf, counts)) {
    stop("Triplet TF-IDF reconstruction changed the ATAC layout.", call. = FALSE)
  }
  object <- .rc_set_assay_matrix(
    object = object,
    assay = atac_assay,
    layer = "data",
    new_data = tfidf
  )
  object <- .rc_align_normalized_assay(object, atac_assay, "ATAC")
  object@misc$regcompass_atac_normalization <- c(list(
    method = "Signac_TFIDF",
    assembly = "single_triplet_reconstruction",
    scope = "cell_type_all_available_cells",
    celltype_col = celltype_col,
    idf_reference = "all cells of the same cell type",
    n_units_by_celltype = vapply(groups, length, integer(1)),
    n_zero_count_peaks_by_celltype = local_zero,
    celltype_local_zero_peak_policy =
      "retain_as_zero_without_passing_to_RunTFIDF",
    tfidf_method = method,
    scale_factor = scale.factor
  ), filtered$diagnostics)
  object
}

.rc_with_layer2_context <- function(
    monitor, parallel, BPPARAM, expression) {
  previous <- list(
    monitor = .rc_layer2_runtime_context$monitor,
    parallel = .rc_layer2_runtime_context$parallel,
    BPPARAM = .rc_layer2_runtime_context$BPPARAM,
    directional_call = .rc_layer2_runtime_context$directional_call
  )
  on.exit({
    .rc_layer2_runtime_context$monitor <- previous$monitor
    .rc_layer2_runtime_context$parallel <- previous$parallel
    .rc_layer2_runtime_context$BPPARAM <- previous$BPPARAM
    .rc_layer2_runtime_context$directional_call <- previous$directional_call
  }, add = TRUE)
  .rc_layer2_runtime_context$monitor <- monitor
  .rc_layer2_runtime_context$parallel <- isTRUE(parallel)
  .rc_layer2_runtime_context$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE
  .rc_layer2_runtime_context$directional_call <- 0L
  force(expression)
}

.rc_run_microcompass_monitored <- function(..., progress_monitor = NULL) {
  args <- list(...)
  parallel <- isTRUE(args$parallel %||% TRUE)
  BPPARAM <- args$BPPARAM %||% NULL
  .rc_step_monitor_event(
    progress_monitor,
    "layer2_engine_start",
    "starting union-model construction and directional scoring",
    context = list(parallel = parallel)
  )
  answer <- .rc_with_layer2_context(
    progress_monitor, parallel, BPPARAM,
    do.call(.rc_original_run_microcompass_hardening, args)
  )
  .rc_step_monitor_event(
    progress_monitor,
    "layer2_engine_complete",
    "completed union-model construction and directional scoring"
  )
  answer
}

.rc_run_microcompass_engine_monitored <- function(
    ..., progress_monitor = NULL) {
  args <- list(...)
  parallel <- isTRUE(args$parallel %||% TRUE)
  BPPARAM <- args$BPPARAM %||% NULL
  .rc_step_monitor_event(
    progress_monitor,
    "layer2_control_start",
    "starting shared-model control scoring"
  )
  answer <- .rc_with_layer2_context(
    progress_monitor, parallel, BPPARAM,
    do.call(.rc_original_run_microcompass_engine_hardening, args)
  )
  .rc_step_monitor_event(
    progress_monitor,
    "layer2_control_complete",
    "completed shared-model control scoring"
  )
  answer
}

.rc_build_medium_specific_union_gem_cache <- function(...) {
  monitor <- .rc_layer2_runtime_context$monitor
  .rc_step_monitor_event(
    monitor,
    "union_gem_start",
    "building or validating medium-specific union GEMs"
  )
  answer <- .rc_original_union_cache_hardening(...)
  .rc_step_monitor_event(
    monitor,
    "union_gem_complete",
    "medium-specific union GEMs are ready",
    context = list(cache_entries = length(answer))
  )
  answer
}

.rc_lp_block_size <- function() {
  value <- suppressWarnings(as.integer(
    getOption("RegCompassR.lp_target_block_size", .rc_lp_target_block_default)
  ))
  if (!is.finite(value) || value < 1L) value <- .rc_lp_target_block_default
  value
}

.rc_with_lp_worker_single_thread <- function(FUN) {
  if (!is.function(FUN)) stop("`FUN` must be a function.", call. = FALSE)
  state <- .rc_set_internal_single_thread()
  old_highs <- Sys.getenv("HIGHS_THREADS", unset = NA_character_)
  Sys.setenv(HIGHS_THREADS = "1")
  on.exit({
    if (is.na(old_highs)) Sys.unsetenv("HIGHS_THREADS") else
      Sys.setenv(HIGHS_THREADS = old_highs)
    .rc_restore_internal_threads(state)
  }, add = TRUE)
  FUN()
}

.rc_directional_feasibility <- function(
    gem, targets, solver = "highs", time_limit = 60,
    flux_threshold = 1e-8) {
  required <- c("reaction_id", "target_direction")
  if (!is.data.frame(targets) || !all(required %in% colnames(targets))) {
    stop("`targets` must contain reaction_id and target_direction.", call. = FALSE)
  }
  if (!nrow(targets)) {
    return(data.frame(
      reaction_id = character(), target_direction = character(),
      feasible = logical(), vmax = numeric(), solver_status = character(),
      stringsAsFactors = FALSE
    ))
  }
  validated <- rc_validate_gem(gem)
  block_size <- .rc_lp_block_size()
  indices <- split(
    seq_len(nrow(targets)),
    ceiling(seq_len(nrow(targets)) / block_size)
  )
  monitor <- .rc_layer2_runtime_context$monitor
  .rc_layer2_runtime_context$directional_call <-
    as.integer(.rc_layer2_runtime_context$directional_call %||% 0L) + 1L
  call_id <- .rc_layer2_runtime_context$directional_call
  .rc_step_monitor_event(
    monitor,
    "directional_feasibility_start",
    "solving directional feasibility in target blocks",
    context = list(call = call_id, targets = nrow(targets), blocks = length(indices))
  )
  run_block <- function(rows) {
    .rc_with_lp_worker_single_thread(function() {
      do.call(rbind, lapply(rows, function(i) {
        reaction <- as.character(targets$reaction_id[[i]])
        direction <- as.character(targets$target_direction[[i]])
        if (!reaction %in% validated$reactions) {
          return(data.frame(
            reaction_id = reaction, target_direction = direction,
            feasible = FALSE, vmax = NA_real_,
            solver_status = "reaction_missing",
            stringsAsFactors = FALSE
          ))
        }
        if (!direction %in% c("forward", "reverse")) {
          return(data.frame(
            reaction_id = reaction, target_direction = direction,
            feasible = FALSE, vmax = 0,
            solver_status = "no_allowed_direction",
            stringsAsFactors = FALSE
          ))
        }
        value <- rc_compass_vmax_directional(
          S = validated$S, lb = validated$lb, ub = validated$ub,
          target_reaction = reaction, direction = direction,
          solver = solver, time_limit = time_limit,
          flux_threshold = flux_threshold
        )
        data.frame(
          reaction_id = reaction,
          target_direction = direction,
          feasible = isTRUE(value$feasible),
          vmax = value$vmax,
          solver_status = value$status,
          stringsAsFactors = FALSE
        )
      }))
    })
  }
  grouped <- rc_parallel_lapply(
    indices,
    run_block,
    BPPARAM = if (isTRUE(.rc_layer2_runtime_context$parallel)) {
      .rc_layer2_runtime_context$BPPARAM
    } else {
      FALSE
    }
  )
  answer <- do.call(rbind, grouped)
  rownames(answer) <- NULL
  .rc_step_monitor_event(
    monitor,
    "directional_feasibility_complete",
    "completed directional feasibility target blocks",
    context = list(call = call_id, targets = nrow(targets), blocks = length(indices))
  )
  answer
}

.rc_build_microcompass_vmax_cache <- function(
    model_cache, mode, model_keys, solver, flux_threshold,
    parallel = TRUE, BPPARAM = NULL) {
  block_size <- .rc_lp_block_size()
  unique_model_keys <- unique(unname(model_keys))
  tasks <- unlist(lapply(unique_model_keys, function(model_key) {
    selected_rows <- names(model_keys)[model_keys == model_key]
    blocks <- split(
      selected_rows,
      ceiling(seq_along(selected_rows) / block_size)
    )
    lapply(blocks, function(row_ids) {
      list(model_key = model_key, row_ids = row_ids)
    })
  }), recursive = FALSE)
  monitor <- .rc_layer2_runtime_context$monitor
  .rc_step_monitor_event(
    monitor,
    "vmax_cache_start",
    "solving structural vmax by model and target block",
    context = list(targets = length(model_cache), blocks = length(tasks))
  )
  grouped <- rc_parallel_lapply(
    tasks,
    function(task) {
      .rc_with_lp_worker_single_thread(function() {
        selected_rows <- as.character(task$row_ids)
        first_entry <- model_cache[[selected_rows[[1L]]]]
        model <- .rc_load_microcompass_model(first_entry, mode)
        values <- lapply(selected_rows, function(row_id) {
          entry <- model_cache[[row_id]]
          rc_compass_vmax_directional(
            S = model$S,
            lb = model$lb,
            ub = model$ub,
            target_reaction = entry$reaction_id,
            direction = entry$target_direction,
            solver = solver,
            flux_threshold = flux_threshold
          )
        })
        names(values) <- selected_rows
        values
      })
    },
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  answer <- .rc_flatten_microcompass_vmax_cache(
    grouped, names(model_cache)
  )
  .rc_step_monitor_event(
    monitor,
    "vmax_cache_complete",
    "completed structural vmax target blocks",
    context = list(targets = length(answer), blocks = length(tasks))
  )
  .rc_step_monitor_event(
    monitor,
    "metacell_scoring_start",
    "starting metacell-specific step-2 LP objectives"
  )
  answer
}

.rc_compact_meta_modules_for_layer2 <- function(meta_modules) {
  embedded <- meta_modules$condition_modules
  if (is.list(embedded) && is.data.frame(embedded$reaction_membership)) {
    meta_modules$condition_modules <- list(
      schema_version = "regcompass_transient_condition_modules_summary_v1",
      embedded = FALSE,
      n_reaction_membership = nrow(embedded$reaction_membership)
    )
    class(meta_modules$condition_modules) <- c(
      "regcompass_external_condition_modules", "list"
    )
  }
  meta_modules
}

rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  meta_modules <- .rc_compact_meta_modules_for_layer2(meta_modules)
  invisible(gc(verbose = FALSE, full = TRUE))
  implementation <- .rc_original_step_layer2_hardening
  evaluation_environment <- new.env(parent = environment(implementation))
  evaluation_environment$rc_run_microcompass <- function(...) {
    monitor <- get("monitor", envir = parent.frame(), inherits = TRUE)
    .rc_run_microcompass_monitored(..., progress_monitor = monitor)
  }
  evaluation_environment$.rc_run_microcompass_engine <- function(...) {
    monitor <- get("monitor", envir = parent.frame(), inherits = TRUE)
    .rc_run_microcompass_engine_monitored(..., progress_monitor = monitor)
  }
  environment(implementation) <- evaluation_environment
  implementation(
    layer1 = layer1,
    meta_modules = meta_modules,
    gem = gem,
    medium_scenarios = medium_scenarios,
    outdir = outdir,
    model_mode = model_mode,
    layer2_args = layer2_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
}

.rc_load_condition_modules <- function(meta_modules) {
  value <- meta_modules$condition_modules
  if (is.list(value) && is.data.frame(value$reaction_membership)) return(value)
  reference <- meta_modules$condition_modules_ref %||% value
  if (!is.list(reference) ||
      !identical(
        reference$schema_version,
        "regcompass_external_condition_modules_v1"
      )) {
    stop("The condition meta-module reference is invalid.", call. = FALSE)
  }
  file <- as.character(reference$file %||% "")
  if (length(file) != 1L || !nzchar(file) || !file.exists(file)) {
    stop("The external condition meta-module file is unavailable.", call. = FALSE)
  }
  checksum <- unname(tools::md5sum(file))
  if (!identical(checksum, as.character(reference$file_checksum))) {
    stop("The external condition meta-module file failed checksum validation.",
         call. = FALSE)
  }
  readRDS(file)
}

rc_regcompass_step_results <- function(
    grn, metacells, meta_modules, layer1, layer2, gem, outdir,
    species = c("auto", "human", "mouse"),
    progress = getOption("RegCompassR.progress", TRUE)) {
  materialized <- meta_modules
  materialized$condition_modules <- .rc_load_condition_modules(meta_modules)
  .rc_original_step_results_hardening(
    grn = grn,
    metacells = metacells,
    meta_modules = materialized,
    layer1 = layer1,
    layer2 = layer2,
    gem = gem,
    outdir = outdir,
    species = species,
    progress = progress
  )
}
