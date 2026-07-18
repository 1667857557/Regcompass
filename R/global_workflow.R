.rc_bind_frames_fill <- function(x) {
  x <- x[vapply(x, function(z) is.data.frame(z) && nrow(z) > 0L, logical(1))]
  if (!length(x)) return(data.frame())
  columns <- unique(unlist(lapply(x, colnames), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(columns, colnames(z))
    for (column in missing) z[[column]] <- NA
    z[, columns, drop = FALSE]
  })
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.rc_cbind_matrix_union <- function(x, fill = NA_real_) {
  x <- x[vapply(x, function(z) !is.null(dim(z)) && ncol(z) > 0L, logical(1))]
  if (!length(x)) return(matrix(numeric(), 0L, 0L))
  rows <- unique(unlist(lapply(x, rownames), use.names = FALSE))
  columns <- unlist(lapply(x, colnames), use.names = FALSE)
  if (anyNA(rows) || any(!nzchar(rows)) || anyNA(columns) || any(!nzchar(columns))) {
    stop("Matrices require non-empty row and column names.", call. = FALSE)
  }
  if (anyDuplicated(columns)) {
    duplicated_ids <- unique(columns[duplicated(columns)])
    stop("Duplicated metacell IDs across strata: ", paste(utils::head(duplicated_ids, 10L), collapse = ", "), call. = FALSE)
  }
  out <- matrix(fill, nrow = length(rows), ncol = length(columns), dimnames = list(rows, columns))
  offset <- 0L
  for (matrix_in in x) {
    column_index <- seq.int(offset + 1L, offset + ncol(matrix_in))
    out[rownames(matrix_in), column_index] <- as.matrix(matrix_in)
    offset <- offset + ncol(matrix_in)
  }
  out
}

.rc_phase_bpparam <- function(workers = NULL, backend = c("auto", "serial", "snow", "multicore")) {
  backend <- match.arg(backend)
  if (identical(backend, "serial")) return(FALSE)
  param <- rc_default_bpparam(workers = workers, backend = backend)
  param %||% FALSE
}

.rc_release_bpparam <- function(param) {
  if (!identical(param, FALSE) && !is.null(param) && requireNamespace("BiocParallel", quietly = TRUE)) {
    try(BiocParallel::bpstop(param), silent = TRUE)
  }
  invisible(gc(verbose = FALSE))
}

.rc_metacell_logcpm <- function(counts, scale_factor = 1e6,
                                 library_size = NULL) {
  counts <- methods::as(counts, "dgCMatrix")
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("`scale_factor` must be one positive finite number.", call. = FALSE)
  }
  normalization_scope <- "input_matrix_library_size"
  if (is.null(library_size)) {
    library_size <- Matrix::colSums(counts)
  } else {
    normalization_scope <- "full_transcriptome_library_size_before_gpr_filter"
    if (!is.null(names(library_size))) {
      missing <- setdiff(colnames(counts), names(library_size))
      if (length(missing)) {
        stop("`library_size` is missing metacells: ",
             paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
      }
      library_size <- library_size[colnames(counts)]
    }
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) ||
      any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop("`library_size` must contain one positive finite value per metacell.",
         call. = FALSE)
  }
  scaled <- counts %*% Matrix::Diagonal(x = scale_factor / library_size)
  answer <- log1p(scaled)
  dimnames(answer) <- dimnames(counts)
  attr(answer, "normalization_scope") <- normalization_scope
  attr(answer, "library_size") <- stats::setNames(
    library_size,
    colnames(counts)
  )
  answer
}


.rc_weighted_quantile <- function(x, weights, probs = 0.5) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- x[keep]
  weights <- weights[keep]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  order_index <- order(x)
  x <- x[order_index]
  weights <- weights[order_index]
  cumulative <- cumsum(weights) / sum(weights)
  vapply(probs, function(probability) {
    if (!is.finite(probability) || probability < 0 || probability > 1) {
      stop("`probs` values must be finite numbers between 0 and 1.", call. = FALSE)
    }
    if (probability <= 0) return(x[[1L]])
    if (probability >= 1) return(x[[length(x)]])
    index <- which(cumulative >= probability)[[1L]]
    x[[index]]
  }, numeric(1))
}

.rc_equal_sample_weights <- function(sample_ids) {
  sample_ids <- trimws(as.character(sample_ids))
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    stop("Sample IDs used for calibration must be non-missing and non-empty.", call. = FALSE)
  }
  counts <- table(sample_ids)
  weights <- 1 / as.numeric(counts[sample_ids])
  names(weights) <- names(sample_ids)
  weights / sum(weights)
}

.rc_weighted_gene_score <- function(
    X, weights, min_scale = 0.05, z_clip = 6,
    mode = c("absolute", "relative"),
    half_saturation = getOption("RegCompassR.cpm_half_saturation", 1)) {
  X <- as.matrix(X)
  if (length(weights) != ncol(X) || any(!is.finite(weights)) ||
      any(weights <= 0)) {
    stop("`weights` must contain one positive finite value per column.",
         call. = FALSE)
  }
  mode <- match.arg(mode)
  if (identical(mode, "absolute")) {
    return(.rc_absolute_activity_score(X, half_saturation))
  }
  centers <- apply(
    X, 1L, .rc_weighted_quantile,
    weights = weights, probs = 0.5
  )
  scales <- vapply(seq_len(nrow(X)), function(i) {
    mad_sigma <- .rc_weighted_quantile(
      abs(X[i, ] - centers[[i]]), weights, probs = 0.5
    ) * 1.4826
    quartiles <- .rc_weighted_quantile(
      X[i, ], weights, probs = c(0.25, 0.75)
    )
    max(mad_sigma, diff(quartiles) / 1.349, min_scale, na.rm = TRUE)
  }, numeric(1))
  z <- sweep(X, 1L, centers, "-")
  z <- sweep(z, 1L, scales, "/")
  z <- pmax(pmin(z, z_clip), -z_clip)
  score <- rc_sigmoid(z)
  finite_range <- apply(X, 1L, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    diff(range(x))
  })
  score[!is.finite(finite_range) | finite_range <= 1e-12, ] <- NA_real_
  dimnames(score) <- dimnames(X)
  attr(score, "score_semantics") <-
    "sample_balanced_within_gene_relative_state"
  score
}


.rc_weighted_q95_calibrate <- function(
    C_raw, weights, eps = 1e-6, n0 = 80,
    unit_meta = NULL, stratum_col = NULL,
    bootstrap = FALSE, B = 500, BPPARAM = NULL) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw))) {
    rownames(C_raw) <- paste0("reaction_", seq_len(nrow(C_raw)))
  }
  if (is.null(colnames(C_raw))) {
    colnames(C_raw) <- paste0("unit_", seq_len(ncol(C_raw)))
  }
  rc_q95_calibrate(
    C_raw = C_raw,
    eps = eps,
    bootstrap = bootstrap,
    B = B,
    BPPARAM = BPPARAM,
    n0 = n0,
    unit_meta = unit_meta,
    stratum_col = stratum_col,
    weights = weights
  )
}


.rc_normalize_calibration_params <- function(
    params, sample_col, condition_col, celltype_col) {
  params <- params %||% list()
  correction <- match.arg(
    params$expression_batch_correction %||% "none",
    c("none", "limma")
  )
  sample_balance <- params$sample_balance %||% TRUE
  q95_bootstrap <- params$q95_bootstrap %||% FALSE
  logical_values <- list(
    sample_balance = sample_balance,
    q95_bootstrap = q95_bootstrap
  )
  invalid_logical <- names(logical_values)[!vapply(
    logical_values,
    function(value) is.logical(value) && length(value) == 1L && !is.na(value),
    logical(1)
  )]
  if (length(invalid_logical)) {
    stop(
      "Calibration switches must be TRUE or FALSE: ",
      paste(invalid_logical, collapse = ", "),
      call. = FALSE
    )
  }
  q95_n0 <- params$q95_n0 %||% 80
  q95_B <- params$q95_B %||% 500L
  if (!is.numeric(q95_n0) || length(q95_n0) != 1L ||
      !is.finite(q95_n0) || q95_n0 < 0) {
    stop("`q95_n0` must be one finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(q95_B) || length(q95_B) != 1L ||
      !is.finite(q95_B) || q95_B < 1) {
    stop("`q95_B` must be one positive finite number.", call. = FALSE)
  }
  q95_stratum_col <- params$q95_stratum_col %||% celltype_col
  if (is.null(q95_stratum_col) || identical(q95_stratum_col, "none")) {
    q95_stratum_col <- NULL
  } else {
    q95_stratum_col <- as.character(q95_stratum_col)
    if (length(q95_stratum_col) != 1L || is.na(q95_stratum_col) ||
        !nzchar(q95_stratum_col)) {
      stop("`q95_stratum_col` must be one metadata column or NULL.",
           call. = FALSE)
    }
  }
  technical_batch_cols <- unique(as.character(
    params$technical_batch_cols %||% character()
  ))
  technical_batch_cols <- technical_batch_cols[
    !is.na(technical_batch_cols) & nzchar(technical_batch_cols)
  ]
  preserve_design_cols <- unique(as.character(
    params$preserve_design_cols %||% c(condition_col, celltype_col)
  ))
  preserve_design_cols <- preserve_design_cols[
    !is.na(preserve_design_cols) & nzchar(preserve_design_cols)
  ]
  list(
    sample_balance = sample_balance,
    sample_balance_col = as.character(
      params$sample_balance_col %||% sample_col
    ),
    expression_batch_correction = correction,
    technical_batch_cols = technical_batch_cols,
    preserve_design_cols = preserve_design_cols,
    q95_n0 = as.numeric(q95_n0),
    q95_stratum_col = q95_stratum_col,
    q95_bootstrap = q95_bootstrap,
    q95_B = as.integer(q95_B)
  )
}


.rc_apply_limma_batch_correction <- function(X, unit_meta, calibration_params,
                                             sample_col) {
  mode <- calibration_params$expression_batch_correction
  if (identical(mode, "none")) {
    return(list(
      matrix = X,
      diagnostics = data.frame(
        method = "none",
        n_batch_levels = NA_integer_,
        batch_columns = "",
        preserve_columns = paste(
          calibration_params$preserve_design_cols,
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    ))
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for expression batch correction.", call. = FALSE)
  }
  batch_cols <- calibration_params$technical_batch_cols
  if (!length(batch_cols)) {
    stop(
      "`technical_batch_cols` is required when `expression_batch_correction = 'limma'`.",
      call. = FALSE
    )
  }
  if (sample_col %in% batch_cols) {
    stop(
      "`sample_id` is a biological replicate and cannot be removed as a technical batch.",
      call. = FALSE
    )
  }
  missing <- setdiff(
    unique(c(batch_cols, calibration_params$preserve_design_cols)),
    colnames(unit_meta)
  )
  if (length(missing)) {
    stop("Batch-correction metadata are missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  batch_frame <- unit_meta[, batch_cols, drop = FALSE]
  if (anyNA(batch_frame)) {
    stop("Technical batch columns contain missing values.", call. = FALSE)
  }
  batch <- do.call(
    interaction,
    c(lapply(batch_frame, as.factor), list(drop = TRUE, lex.order = TRUE))
  )
  n_batch_levels <- nlevels(batch)
  if (n_batch_levels < 2L) {
    warning(
      "Only one technical batch level was detected; limma correction was skipped.",
      call. = FALSE
    )
    return(list(
      matrix = X,
      diagnostics = data.frame(
        method = "limma_skipped_single_level",
        n_batch_levels = n_batch_levels,
        batch_columns = paste(batch_cols, collapse = ";"),
        preserve_columns = paste(
          calibration_params$preserve_design_cols,
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    ))
  }
  preserve_cols <- calibration_params$preserve_design_cols
  preserve_cols <- preserve_cols[vapply(
    unit_meta[, preserve_cols, drop = FALSE],
    function(value) length(unique(value[!is.na(value)])) > 1L,
    logical(1)
  )]
  design <- if (length(preserve_cols)) {
    stats::model.matrix(stats::reformulate(preserve_cols), data = unit_meta)
  } else {
    matrix(1, nrow = nrow(unit_meta), ncol = 1L)
  }
  confounding_data <- unit_meta
  confounding_data$.rc_technical_batch <- batch
  full_terms <- c(preserve_cols, ".rc_technical_batch")
  full_design <- stats::model.matrix(
    stats::reformulate(full_terms),
    data = confounding_data
  )
  if (qr(full_design)$rank < ncol(full_design)) {
    stop(
      "Technical batch is confounded with the preserved biological design; correction is not identifiable.",
      call. = FALSE
    )
  }
  corrected <- limma::removeBatchEffect(
    as.matrix(X),
    batch = batch,
    design = design
  )
  dimnames(corrected) <- dimnames(X)
  list(
    matrix = corrected,
    diagnostics = data.frame(
      method = "limma_removeBatchEffect",
      n_batch_levels = n_batch_levels,
      batch_columns = paste(batch_cols, collapse = ";"),
      preserve_columns = paste(preserve_cols, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

.rc_local_fastcore_rows <- function(model, membership, core_ids) {
  reactions <- colnames(model$S)
  template <- membership[rep(1L, length(reactions)), , drop = FALSE]
  template$reaction_id <- reactions
  matched <- match(reactions, as.character(membership$reaction_id))
  existing <- !is.na(matched)
  if (any(existing)) {
    template[existing, colnames(membership)] <- membership[
      matched[existing],
      colnames(membership),
      drop = FALSE
    ]
  }
  support <- character()
  if (!is.null(model$reaction_meta) &&
      all(c("reaction_id", "fastcore_support") %in% colnames(model$reaction_meta))) {
    support <- as.character(model$reaction_meta$reaction_id[
      model$reaction_meta$fastcore_support %in% TRUE
    ])
  }
  template$is_core <- reactions %in% core_ids
  template$biological_meta_module_member <- reactions %in%
    as.character(membership$reaction_id)
  template$local_fastcore_support <- reactions %in% support
  previous_stage <- if ("inclusion_stage" %in% colnames(template)) {
    as.character(template$inclusion_stage)
  } else {
    rep(NA_character_, nrow(template))
  }
  previous_stage[is.na(previous_stage) | !nzchar(previous_stage)] <-
    "biological_meta_module_member"
  template$inclusion_stage <- ifelse(
    template$local_fastcore_support,
    "local_fastcore_support",
    previous_stage
  )
  template
}

.rc_complete_stratum_meta_modules <- function(meta_modules, gem, outdir,
                                              local_fastcore_args = list()) {
  defaults <- list(
    enabled = TRUE,
    target_direction = "both",
    solver = "highs",
    time_limit = 300,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE,
    save_models = TRUE
  )
  defaults[names(local_fastcore_args)] <- NULL
  args <- c(defaults, local_fastcore_args)
  if (!isTRUE(args$enabled)) {
    return(list(
      completed_reaction_membership = meta_modules$reaction_membership,
      summary = data.frame(),
      diagnostics = data.frame(),
      completion_iterations = data.frame(),
      parent_scope = "disabled"
    ))
  }
  membership <- meta_modules$reaction_membership
  core <- meta_modules$core_gene_reaction
  required <- c("sample_id", "module_id", "reaction_id")
  if (!is.data.frame(membership) || !all(required %in% colnames(membership))) {
    stop("Meta-module reaction membership is incomplete before local FASTCORE.", call. = FALSE)
  }
  if (!is.data.frame(core) || !all(required %in% colnames(core))) {
    stop("Meta-module core reactions are incomplete before local FASTCORE.", call. = FALSE)
  }
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  parent <- .rc_fastcore_parent(
    gem,
    medium_table = NULL,
    condition = NULL,
    solver = args$solver,
    time_limit = args$time_limit,
    fastcore_epsilon = args$fastcore_epsilon
  )
  module_keys <- unique(membership[, c("sample_id", "module_id"), drop = FALSE])
  completed_rows <- vector("list", nrow(module_keys))
  summaries <- vector("list", nrow(module_keys))
  diagnostics <- vector("list", nrow(module_keys))
  iterations <- vector("list", nrow(module_keys))
  model_dir <- file.path(outdir, "local_fastcore_models")
  if (isTRUE(args$save_models)) {
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  }
  for (index in seq_len(nrow(module_keys))) {
    sample_id <- as.character(module_keys$sample_id[[index]])
    module_id <- as.character(module_keys$module_id[[index]])
    in_module <- as.character(membership$sample_id) == sample_id &
      as.character(membership$module_id) == module_id
    membership_i <- membership[in_module, , drop = FALSE]
    in_core <- as.character(core$sample_id) == sample_id &
      as.character(core$module_id) == module_id
    core_i <- core[in_core, , drop = FALSE]
    if (!nrow(core_i)) {
      stop("No core reactions remain for local meta-module `", module_id, "`.",
           call. = FALSE)
    }
    model <- .rc_complete_meta_module(
      gem = gem,
      reaction_membership = membership_i,
      core_reactions = core_i,
      sample_id = sample_id,
      module_id = module_id,
      medium_table = NULL,
      condition = NULL,
      parent_gem = parent,
      target_direction = args$target_direction,
      solver = args$solver,
      time_limit = args$time_limit,
      fastcore_epsilon = args$fastcore_epsilon,
      max_support_reactions = args$max_support_reactions,
      strict = args$strict
    )
    completed_rows[[index]] <- .rc_local_fastcore_rows(
      model,
      membership_i,
      unique(as.character(core_i$reaction_id))
    )
    support_ids <- model$reaction_meta$reaction_id[
      model$reaction_meta$fastcore_support %in% TRUE
    ]
    summaries[[index]] <- data.frame(
      sample_id = sample_id,
      module_id = module_id,
      n_biological_reactions = nrow(membership_i),
      n_core_reactions = length(unique(as.character(core_i$reaction_id))),
      n_local_fastcore_support = length(unique(as.character(support_ids))),
      n_completed_reactions = ncol(model$S),
      target_status = model$target_status,
      parent_scope = "unconstrained_shared_parent_with_fastcc",
      stringsAsFactors = FALSE
    )
    if (nrow(model$closure_diagnostics)) {
      diagnostic_i <- model$closure_diagnostics
      diagnostic_i$sample_id <- sample_id
      diagnostic_i$module_id <- module_id
      diagnostics[[index]] <- diagnostic_i
    }
    if (nrow(model$completion_iterations)) {
      iteration_i <- model$completion_iterations
      iteration_i$sample_id <- sample_id
      iteration_i$module_id <- module_id
      iterations[[index]] <- iteration_i
    }
    if (isTRUE(args$save_models)) {
      safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", module_id)
      saveRDS(model, file.path(model_dir, paste0(safe_id, ".rds")))
    }
  }
  list(
    completed_reaction_membership = .rc_bind_frames_fill(completed_rows),
    summary = .rc_bind_frames_fill(summaries),
    diagnostics = .rc_bind_frames_fill(diagnostics),
    completion_iterations = .rc_bind_frames_fill(iterations),
    parent_scope = "unconstrained_shared_parent_with_fastcc"
  )
}

.rc_run_regcompass_stratum <- function(object, group_id, group_cols, gem, outdir,
                                        pfm, genome, fragment_files = NULL,
                                        sample_col = "sample_id",
                                        condition_col = "condition",
                                        celltype_col = "cell_type",
                                        rna_assay = "RNA",
                                        atac_assay = "ATAC",
                                        metacell_args = list(),
                                        layer1_args = list(),
                                        pando_args = list()) {
  allowed_layer1_args <- c(
    "promiscuity_mode", "and_method", "or_method", "tau",
    "local_fastcore", "local_fastcore_args",
    "sample_balance", "sample_balance_col",
    "expression_batch_correction", "technical_batch_cols",
    "preserve_design_cols", "q95_n0", "q95_stratum_col",
    "q95_bootstrap", "q95_B"
  )
  unsupported_layer1_args <- setdiff(names(layer1_args), allowed_layer1_args)
  if (length(unsupported_layer1_args)) {
    stop("Unsupported `layer1_args`: ", paste(unsupported_layer1_args, collapse = ", "), call. = FALSE)
  }
  meta <- object@meta.data
  ids <- rc_make_stratum_id(meta, group_cols)
  cells <- rownames(meta)[ids == group_id]
  if (!length(cells)) {
    stop("No cells found for stratum: ", group_id, call. = FALSE)
  }
  one <- subset(object, cells = cells)
  stratum_dir <- file.path(
    outdir,
    gsub("[^A-Za-z0-9_.-]+", "_", group_id)
  )
  dir.create(stratum_dir, recursive = TRUE, showWarnings = FALSE)
  or_method <- match.arg(
    layer1_args$or_method %||% "max",
    c("max", "sum_sqrtK", "prob_or", "sum")
  )
  capacity_params <- list(
    promiscuity_mode = match.arg(
      layer1_args$promiscuity_mode %||% "none",
      c("none", "sqrt", "linear")
    ),
    and_method = match.arg(
      layer1_args$and_method %||% "min",
      c("min", "boltzmann", "mean")
    ),
    tau = layer1_args$tau %||% 0.20,
    or_method = or_method
  )
  fragment_aggregation_enabled <- !identical(fragment_files, FALSE)
  metacell_defaults <- list(
    object = one,
    outdir = file.path(stratum_dir, "01_metacells"),
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = if (isTRUE(fragment_aggregation_enabled)) fragment_files else FALSE,
    save_metacell_object = TRUE,
    save_counts = TRUE,
    save_fragments = fragment_aggregation_enabled,
    require_fragment_aggregation = fragment_aggregation_enabled,
    fragment_aggregation_backend = if (isTRUE(fragment_aggregation_enabled)) "regcompass" else "none",
    BPPARAM = FALSE,
    on_stratum_error = "stop"
  )
  reserved <- intersect(names(metacell_args), names(metacell_defaults))
  reserved <- setdiff(
    reserved,
    c(
      "rna_reduction", "atac_reduction", "rna_dims", "atac_dims",
      "gamma", "seed", "min_cells_per_stratum",
      "min_metacell_size", "min_metacells_per_stratum",
      "label_col", "bgzip_path", "tabix_path",
      "fragment_nb_cl", "overwrite"
    )
  )
  if (length(reserved)) {
    stop(
      "`metacell_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }
  metacell_defaults[names(metacell_args)] <- NULL
  metacells <- do.call(
    rc_make_supercell2_metacells,
    c(metacell_defaults, metacell_args)
  )
  minimum_metacells <- as.integer(pando_args$min_metacells %||% 20L)
  if (nrow(metacells$metacell_meta) < minimum_metacells) {
    return(list(
      group_id = group_id,
      status = "skipped_too_few_metacells",
      artifact_file = NA_character_,
      n_cells = length(cells),
      n_metacells = nrow(metacells$metacell_meta),
      error_class = NA_character_,
      error_message = paste0("Produced fewer than ", minimum_metacells, " metacells.")
    ))
  }
  metacell_object <- rc_load_or_merge_metacell_objects(
    metacells$metacell_objects,
    fragment_manifest = metacells$fragment_manifest,
    metacell_meta = metacells$metacell_meta,
    fragment_files = metacells$fragment_files,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    require_complete_fragments = fragment_aggregation_enabled
  )
  pando_args$BPPARAM <- NULL
  pando_args$save_sample_metacell_objects <- NULL
  pando_infer_args <- pando_args$pando_infer_args %||% list()
  pando_infer_args$parallel <- FALSE
  pando_args$pando_infer_args <- pando_infer_args
  pando_args$group_cols <- NULL
  pando_args$sample_col <- NULL
  pando_args$condition_col <- NULL
  pando_args$celltype_col <- NULL
  pando_outdir <- file.path(stratum_dir, "02_pando_meta_modules")
  pando_defaults <- list(
    metacell_object = metacell_object,
    gem = gem,
    outdir = pando_outdir,
    pfm = pfm,
    genome = genome,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    group_cols = group_cols,
    single_cell_genes = rownames(.rc_get_assay_counts(one, rna_assay)),
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    min_metacells = minimum_metacells,
    save_sample_metacell_objects = TRUE,
    BPPARAM = FALSE,
    on_sample_error = "stop"
  )
  pando_defaults[names(pando_args)] <- NULL
  meta_modules <- do.call(
    rc_run_pando_meta_modules,
    c(pando_defaults, pando_args)
  )
  pando_object_file <- file.path(
    pando_outdir,
    "sample_metacell_objects",
    paste0(gsub("[^A-Za-z0-9_.-]+", "_", group_id), ".rds")
  )
  if (!file.exists(pando_object_file)) {
    stop(
      "Pando did not save its normalized metacell object: ",
      pando_object_file,
      call. = FALSE
    )
  }
  pando_object <- readRDS(pando_object_file)
  pando_confidence <- .rc_pando_reaction_confidence(
    meta_modules,
    pando_object,
    gem,
    atac_assay = atac_assay,
    rna_assay = rna_assay
  )
  saveRDS(
    pando_confidence$gene_confidence,
    file.path(pando_outdir, "pando_gene_confidence.rds")
  )
  saveRDS(
    pando_confidence$reaction_confidence_matrix,
    file.path(pando_outdir, "pando_reaction_confidence_matrix.rds")
  )
  .rc_mm_write_tsv_gz(
    pando_confidence$gene_confidence_diagnostics,
    file.path(pando_outdir, "pando_gene_confidence_diagnostics.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    pando_confidence$reaction_confidence,
    file.path(pando_outdir, "pando_reaction_confidence.tsv.gz")
  )
  meta_modules$gene_confidence <- pando_confidence$gene_confidence
  meta_modules$gene_confidence_diagnostics <-
    pando_confidence$gene_confidence_diagnostics
  meta_modules$reaction_confidence <-
    pando_confidence$reaction_confidence
  meta_modules$reaction_confidence_matrix <-
    pando_confidence$reaction_confidence_matrix
  meta_modules$confidence_source <- pando_confidence$confidence_source
  meta_modules$biological_reaction_membership <- meta_modules$reaction_membership

  local_fastcore_args <- layer1_args$local_fastcore_args %||% list()
  local_fastcore_args$enabled <- layer1_args$local_fastcore %||%
    local_fastcore_args$enabled %||% TRUE
  local_completion <- .rc_complete_stratum_meta_modules(
    meta_modules,
    gem,
    outdir = file.path(stratum_dir, "03_local_fastcore"),
    local_fastcore_args = local_fastcore_args
  )
  meta_modules$local_completed_reaction_membership <-
    local_completion$completed_reaction_membership
  meta_modules$local_fastcore_summary <- local_completion$summary
  meta_modules$local_fastcore_diagnostics <- local_completion$diagnostics
  meta_modules$local_fastcore_completion_iterations <-
    local_completion$completion_iterations
  meta_modules$local_fastcore_parent_scope <- local_completion$parent_scope

  gpr_genes <- toupper(rc_metabolic_gpr_genes(gem$gpr_table))
  full_library_size <- Matrix::colSums(metacells$rna_counts)
  if (any(!is.finite(full_library_size)) || any(full_library_size <= 0)) {
    stop("Every metacell must have a positive finite full RNA library size.",
         call. = FALSE)
  }
  rna_counts <- metacells$rna_counts[
    toupper(rownames(metacells$rna_counts)) %in% gpr_genes,
    ,
    drop = FALSE
  ]
  if (!nrow(rna_counts)) {
    stop(
      "No Human-GEM GPR genes were retained in metacell RNA counts.",
      call. = FALSE
    )
  }
  rna_logcpm <- .rc_metacell_logcpm(
    rna_counts,
    library_size = full_library_size[colnames(rna_counts)]
  )
  unit_meta <- metacells$metacell_meta
  id_col <- if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else {
    NULL
  }
  if (is.null(id_col)) {
    stop("Metacell metadata lacks metacell_id/pool_id.", call. = FALSE)
  }
  if (!"pool_id" %in% colnames(unit_meta)) {
    unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  }
  if (!"unit_id" %in% colnames(unit_meta)) {
    unit_meta$unit_id <- unit_meta$pool_id
  }
  unit_meta <- unit_meta[
    match(colnames(rna_logcpm), as.character(unit_meta$pool_id)),
    ,
    drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Metacell metadata are incomplete after Pando.", call. = FALSE)
  }
  reaction_confidence <- as.matrix(
    pando_confidence$reaction_confidence_matrix
  )
  missing_confidence_units <- setdiff(
    colnames(rna_logcpm),
    colnames(reaction_confidence)
  )
  if (length(missing_confidence_units)) {
    stop(
      "Pando reaction confidence is missing metacells: ",
      paste(utils::head(missing_confidence_units, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  reaction_confidence <- reaction_confidence[
    ,
    colnames(rna_logcpm),
    drop = FALSE
  ]
  layer1 <- list(
    schema_version = "regcompass_stratum_evidence_v2",
    rna_metacell_logcpm = rna_logcpm,
    reaction_confidence = reaction_confidence,
    reaction_confidence_source =
      "pando_internal_peak_gene_accessibility",
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    strict_group_id = group_id
  )
  calibration_params <- .rc_normalize_calibration_params(
    layer1_args,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  artifact <- list(
    schema_version = "regcompass_stratum_v3",
    group_id = group_id,
    group_cols = group_cols,
    capacity_params = capacity_params,
    calibration_params = calibration_params,
    layer1 = layer1,
    grn_meta_modules = meta_modules,
    metacell_meta = unit_meta,
    metacell_dir = stratum_dir
  )
  artifact_file <- file.path(stratum_dir, "stratum_result.rds")
  saveRDS(artifact, artifact_file)
  list(
    group_id = group_id,
    status = "ok",
    artifact_file = artifact_file,
    n_cells = length(cells),
    n_metacells = nrow(unit_meta),
    error_class = NA_character_,
    error_message = NA_character_
  )
}

.rc_validate_stratum_artifact_contract <- function(artifact) {
  if (!identical(artifact$schema_version, "regcompass_stratum_v3")) {
    return("schema_version")
  }
  if (!is.list(artifact$layer1) || !is.list(artifact$grn_meta_modules)) {
    return("layer1_or_grn_missing")
  }
  layer1 <- artifact$layer1
  if (is.null(layer1$rna_metacell_logcpm) ||
      is.null(layer1$reaction_confidence) ||
      !is.data.frame(layer1$unit_meta)) {
    return("layer1_outputs_missing")
  }
  rna <- as.matrix(layer1$rna_metacell_logcpm)
  confidence <- as.matrix(layer1$reaction_confidence)
  if (!nrow(rna) || !ncol(rna)) return("empty_rna_metacell_logcpm")
  if (!nrow(confidence) || !ncol(confidence)) return("empty_reaction_confidence")
  unit_meta <- layer1$unit_meta
  id_col <- if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else {
    return("unit_meta_id_missing")
  }
  if (!identical(colnames(rna), as.character(unit_meta[[id_col]]))) {
    return("rna_unit_meta_mismatch")
  }
  missing_confidence <- setdiff(colnames(rna), colnames(confidence))
  if (length(missing_confidence)) return("reaction_confidence_unit_mismatch")
  grn <- artifact$grn_meta_modules
  required_grn <- c("core_gene_reaction", "reaction_membership")
  missing_grn <- required_grn[!vapply(grn[required_grn], is.data.frame, logical(1))]
  if (length(missing_grn)) return(paste0("grn_missing_", paste(missing_grn, collapse = "_")))
  if (!"reaction_id" %in% colnames(grn$core_gene_reaction) ||
      !"reaction_id" %in% colnames(grn$reaction_membership)) {
    return("grn_reaction_id_missing")
  }
  "ok"
}

.rc_merge_stratum_meta_modules <- function(artifacts) {
  names_to_merge <- c(
    "sample_status", "tf_peak_gene_all", "tf_peak_gene_significant",
    "metabolic_gene_nodes", "metabolic_gene_edges", "core_gene_reaction",
    "reaction_membership", "meta_module_summary"
  )
  out <- lapply(names_to_merge, function(name) {
    .rc_bind_frames_fill(lapply(
      artifacts,
      function(artifact) artifact$grn_meta_modules[[name]]
    ))
  })
  names(out) <- names_to_merge

  core <- out$core_gene_reaction
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  core_ids <- unique(as.character(core$reaction_id))
  core_ids <- core_ids[!is.na(core_ids) & nzchar(core_ids)]

  completed <- .rc_bind_frames_fill(lapply(artifacts, function(artifact) {
    artifact$grn_meta_modules$local_completed_reaction_membership %||%
      artifact$grn_meta_modules$reaction_membership
  }))
  if (!length(core_ids) || !nrow(completed)) {
    stop("No completed global meta-module reactions were produced.", call. = FALSE)
  }

  biological <- out$reaction_membership
  biological_ids <- unique(as.character(biological$reaction_id))
  completed_ids <- unique(as.character(completed$reaction_id))
  completed_ids <- completed_ids[!is.na(completed_ids) & nzchar(completed_ids)]

  out$biological_reaction_membership <- biological
  out$local_completed_reaction_membership <- completed
  out$local_fastcore_summary <- .rc_bind_frames_fill(lapply(
    artifacts,
    function(artifact) artifact$grn_meta_modules$local_fastcore_summary
  ))
  out$local_fastcore_diagnostics <- .rc_bind_frames_fill(lapply(
    artifacts,
    function(artifact) artifact$grn_meta_modules$local_fastcore_diagnostics
  ))
  out$local_fastcore_completion_iterations <- .rc_bind_frames_fill(lapply(
    artifacts,
    function(artifact) artifact$grn_meta_modules$local_fastcore_completion_iterations
  ))
  out$global_core_reactions <- data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    reaction_id = core_ids,
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  out$global_reaction_membership <- data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    reaction_id = completed_ids,
    is_core = completed_ids %in% core_ids,
    inclusion_stage = ifelse(
      completed_ids %in% core_ids,
      "global_union_core",
      ifelse(
        completed_ids %in% biological_ids,
        "global_union_biological_member",
        "global_union_local_fastcore_support"
      )
    ),
    stringsAsFactors = FALSE
  )
  out$schema_version <- "regcompass_global_meta_module_v3"
  out$source_group_ids <- unique(vapply(
    artifacts,
    function(artifact) as.character(artifact$group_id),
    character(1)
  ))
  out$global_union_source <-
    "deduplicated_local_fastcore_completed_meta_modules"
  out
}

.rc_merge_stratum_layer1 <- function(artifacts, gem, single_cell_genes,
                                      sample_col, condition_col, celltype_col) {
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  expression_list <- lapply(
    artifacts,
    function(artifact) artifact$layer1$rna_metacell_logcpm
  )
  if (any(vapply(expression_list, is.null, logical(1)))) {
    stop(
      "Every upstream artifact must contain metacell RNA logCPM for global capacity recomputation.",
      call. = FALSE
    )
  }
  expression_list <- lapply(expression_list, function(expression) {
    keep <- tolower(rownames(expression)) %in% gpr_genes
    as.matrix(expression[keep, , drop = FALSE])
  })
  common_genes <- Reduce(intersect, lapply(expression_list, rownames))
  common_genes <- rownames(expression_list[[1L]])[
    rownames(expression_list[[1L]]) %in% common_genes
  ]
  if (!length(common_genes)) {
    stop("No common GPR genes remain across upstream metacell artifacts.", call. = FALSE)
  }
  expression_list <- lapply(
    expression_list,
    function(expression) expression[common_genes, , drop = FALSE]
  )
  rna_logcpm_uncorrected <- .rc_cbind_matrix_union(expression_list)
  unit_meta <- .rc_bind_frames_fill(lapply(
    artifacts,
    function(artifact) artifact$layer1$unit_meta
  ))
  id_col <- if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else {
    NULL
  }
  if (is.null(id_col)) {
    stop("Merged metacell metadata lacks pool_id/metacell_id.", call. = FALSE)
  }
  unit_meta <- unit_meta[
    match(colnames(rna_logcpm_uncorrected), as.character(unit_meta[[id_col]])),
    ,
    drop = FALSE
  ]
  if (anyNA(unit_meta[[id_col]])) {
    stop("Merged metacell metadata are incomplete.", call. = FALSE)
  }
  if (!"pool_id" %in% colnames(unit_meta)) {
    unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  }
  if (!"unit_id" %in% colnames(unit_meta)) unit_meta$unit_id <- unit_meta$pool_id
  unit_meta$stratum_id <- rc_make_stratum_id(
    unit_meta,
    c(condition_col, sample_col, celltype_col)
  )
  parameter_keys <- vapply(artifacts, function(artifact) {
    params <- artifact$capacity_params
    if (is.null(params)) return(NA_character_)
    paste(
      params$promiscuity_mode,
      params$and_method,
      params$tau,
      params$or_method,
      sep = "|"
    )
  }, character(1))
  if (anyNA(parameter_keys) || length(unique(parameter_keys)) != 1L) {
    stop(
      "Upstream artifacts must use one identical Layer 1 capacity parameter set.",
      call. = FALSE
    )
  }
  capacity_params <- artifacts[[1L]]$capacity_params
  calibration_list <- lapply(artifacts, function(artifact) {
    .rc_normalize_calibration_params(
      artifact$calibration_params,
      sample_col = sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col
    )
  })
  calibration_params <- calibration_list[[1L]]
  identical_calibration <- vapply(
    calibration_list,
    identical,
    logical(1),
    y = calibration_params
  )
  if (!all(identical_calibration)) {
    stop("Upstream artifacts use inconsistent global calibration settings.", call. = FALSE)
  }
  correction <- .rc_apply_limma_batch_correction(
    rna_logcpm_uncorrected,
    unit_meta,
    calibration_params,
    sample_col = sample_col
  )
  rna_logcpm <- correction$matrix
  balance_col <- calibration_params$sample_balance_col
  if (!balance_col %in% colnames(unit_meta)) {
    stop("Sample-balance metadata are missing column: ", balance_col, call. = FALSE)
  }
  weights <- if (isTRUE(calibration_params$sample_balance)) {
    .rc_equal_sample_weights(unit_meta[[balance_col]])
  } else {
    rep(1 / ncol(rna_logcpm), ncol(rna_logcpm))
  }
  names(weights) <- colnames(rna_logcpm)
  global_gene_score <- if (isTRUE(calibration_params$sample_balance)) {
    .rc_weighted_gene_score(rna_logcpm, weights)
  } else {
    rc_gene_score(rna_logcpm)
  }
  C_raw <- rc_reaction_capacity(
    parsed,
    global_gene_score,
    promiscuity_mode = capacity_params$promiscuity_mode,
    tau = capacity_params$tau,
    and_method = capacity_params$and_method,
    or_method = capacity_params$or_method,
    BPPARAM = FALSE
  )
  q95_stratum_col <- calibration_params$q95_stratum_col
  if (!is.null(q95_stratum_col) && !q95_stratum_col %in% colnames(unit_meta)) {
    stop("Q95 stratum metadata are missing column: ", q95_stratum_col,
         call. = FALSE)
  }
  calibrated <- rc_q95_calibrate(
    C_raw = C_raw,
    bootstrap = calibration_params$q95_bootstrap,
    B = calibration_params$q95_B,
    BPPARAM = FALSE,
    n0 = calibration_params$q95_n0,
    unit_meta = unit_meta,
    stratum_col = q95_stratum_col,
    weights = if (isTRUE(calibration_params$sample_balance)) weights else NULL
  )
  confidence_list <- lapply(artifacts, function(artifact) {
    value <- artifact$layer1$reaction_confidence
    if (is.null(value)) {
      stop(
        "Every upstream artifact must contain Pando-derived reaction confidence.",
        call. = FALSE
      )
    }
    as.matrix(value)
  })
  reaction_confidence <- .rc_cbind_matrix_union(confidence_list)
  reaction_confidence <- rc_align_layer2_evidence(
    reaction_confidence,
    rownames(C_raw),
    NA_real_
  )
  missing_confidence_units <- setdiff(
    colnames(C_raw),
    colnames(reaction_confidence)
  )
  if (length(missing_confidence_units)) {
    stop(
      "Reaction confidence is missing metacells after global alignment.",
      call. = FALSE
    )
  }
  reaction_confidence <- reaction_confidence[
    ,
    colnames(C_raw),
    drop = FALSE
  ]
  calibration_scope <- paste(
    if (isTRUE(calibration_params$sample_balance)) {
      "equal_sample_weighted"
    } else {
      "equal_metacell_weighted"
    },
    calibration_params$expression_batch_correction,
    "global_gene_score_and_reaction_q95",
    sep = "_"
  )
  list(
    schema_version = "regcompass_global_layer1_v3",
    C_or_raw = C_raw,
    C_raw = C_raw,
    reaction_capacity_L1 = C_raw,
    C_rel = calibrated$C_rel,
    C_abs = calibrated$C_abs,
    C_within_reaction_relative = calibrated$C_within_reaction_relative,
    reaction_confidence = reaction_confidence,
    q95_diagnostics = calibrated$Q,
    capacity_calibration_scope = calibration_scope,
    reaction_confidence_source = "pando_signed_tf_peak_gene_regulatory_support",
    capacity_params = capacity_params,
    calibration_params = calibration_params,
    sample_balance_weights = weights,
    expression_batch_diagnostics = correction$diagnostics,
    rna_metacell_logcpm_uncorrected = rna_logcpm_uncorrected,
    rna_metacell_logcpm = rna_logcpm,
    global_gene_score = global_gene_score,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, tolower(single_cell_genes)),
    parsed_gpr = parsed,
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "metacell",
    strict_stratum_cols = c(condition_col, sample_col, celltype_col)
  )
}
