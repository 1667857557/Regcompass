# Route Stage 1 Pando arguments and schedule condition/standard Pando work.

.rc_stage1_pando_working_object <- function(
    object, rna_assay = "RNA", atac_assay = "ATAC") {
  if (!inherits(object, "Seurat")) {
    stop("Stage 1 Pando working data must inherit from Seurat.", call. = FALSE)
  }
  assays <- unique(c(as.character(rna_assay), as.character(atac_assay)))
  missing <- setdiff(assays, .rc_seurat_assay_names(object))
  if (length(missing)) {
    stop(
      "Stage 1 Pando working data are missing assays: ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  SeuratObject::DefaultAssay(object) <- rna_assay
  slim <- Seurat::DietSeurat(
    object = object, counts = TRUE, data = TRUE, scale.data = FALSE,
    assays = assays, dimreducs = NULL, graphs = NULL, misc = FALSE
  )
  SeuratObject::DefaultAssay(slim) <- rna_assay
  slim
}

.rc_stage1_metabolic_detection_threshold <- 0.20
.rc_stage1_tf_detection_threshold <- 0.05
.rc_stage1_peak_detection_threshold <- 0.05

.rc_stage1_condition_detection <- function(
    counts, metadata, condition_col, features = rownames(counts),
    keep_matrix = FALSE) {
  if (is.null(dim(counts)) || length(dim(counts)) != 2L ||
      is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Stage 1 detection requires a named feature-by-cell count matrix.",
         call. = FALSE)
  }
  if (!is.data.frame(metadata) || is.null(rownames(metadata))) {
    stop("Stage 1 detection requires row-named cell metadata.", call. = FALSE)
  }
  if (!is.character(condition_col) || length(condition_col) != 1L ||
      is.na(condition_col) || !nzchar(trimws(condition_col)) ||
      !condition_col %in% colnames(metadata)) {
    stop("Stage 1 detection requires a valid condition column.", call. = FALSE)
  }
  features <- unique(as.character(features))
  features <- features[!is.na(features) & nzchar(features)]
  missing_features <- setdiff(features, rownames(counts))
  if (length(missing_features)) {
    stop("Stage 1 detection features are absent from the count matrix: ",
         paste(utils::head(missing_features, 10L), collapse = ", "),
         call. = FALSE)
  }
  cells <- as.character(colnames(counts))
  missing_cells <- setdiff(cells, rownames(metadata))
  if (length(missing_cells)) {
    stop("Stage 1 detection metadata are missing count-matrix cells.",
         call. = FALSE)
  }
  meta <- metadata[match(cells, rownames(metadata)), , drop = FALSE]
  condition <- trimws(as.character(meta[[condition_col]]))
  invalid <- is.na(meta[[condition_col]]) | !nzchar(condition)
  if (any(invalid)) {
    stop("Stage 1 detection found missing condition labels.", call. = FALSE)
  }
  condition_levels <- unique(condition)
  max_detection <- stats::setNames(rep(0, length(features)), features)
  condition_detection <- if (isTRUE(keep_matrix)) {
    matrix(
      0, nrow = length(features), ncol = length(condition_levels),
      dimnames = list(features, condition_levels)
    )
  } else NULL
  for (i in seq_along(condition_levels)) {
    level <- condition_levels[[i]]
    cells_use <- cells[condition == level]
    detected <- Matrix::rowSums(
      counts[features, cells_use, drop = FALSE] > 0
    )
    fraction <- as.numeric(detected) / length(cells_use)
    max_detection <- pmax(max_detection, fraction)
    if (isTRUE(keep_matrix)) condition_detection[, i] <- fraction
  }
  list(
    max_detection = max_detection,
    condition_detection = condition_detection,
    condition_levels = condition_levels
  )
}

.rc_stage1_filter_gem_metabolic_genes <- function(
    object, gem, condition_col, rna_assay = "RNA", cell_type = NULL,
    threshold = .rc_stage1_metabolic_detection_threshold) {
  if (!inherits(object, "Seurat")) {
    stop("Stage 1 metabolic target filtering requires a Seurat object.",
         call. = FALSE)
  }
  if (!is.list(gem)) {
    stop("Stage 1 metabolic target filtering requires a GEM list.",
         call. = FALSE)
  }
  threshold <- suppressWarnings(as.numeric(threshold))
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold < 0 || threshold > 1) {
    stop("Stage 1 metabolic detection threshold must be in [0, 1].",
         call. = FALSE)
  }

  counts <- .rc_get_assay_counts(object, rna_assay)
  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(counts)
  target_upper <- intersect(toupper(rna_genes), toupper(metabolic_genes))
  candidate_genes <- rna_genes[toupper(rna_genes) %in% target_upper]
  if (!length(candidate_genes)) {
    stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)
  }

  detection <- .rc_stage1_condition_detection(
    counts = counts, metadata = object@meta.data,
    condition_col = condition_col, features = candidate_genes,
    keep_matrix = TRUE
  )
  retained <- candidate_genes[detection$max_detection >= threshold]
  label <- if (is.null(cell_type) || !length(cell_type) ||
      is.na(cell_type[[1L]]) || !nzchar(trimws(as.character(cell_type[[1L]])))) {
    "<unspecified>"
  } else {
    trimws(as.character(cell_type[[1L]]))
  }
  if (!length(retained)) {
    stop(
      "No GEM metabolic gene reaches RNA detection >= ",
      format(100 * threshold, trim = TRUE),
      "% (counts > 0) in any retained condition for cell type `",
      label, "`.", call. = FALSE
    )
  }

  filtered <- gem
  filtered$metabolic_genes <- retained
  attr(filtered, "regcompass_stage1_metabolic_detection") <- list(
    threshold = threshold,
    positive_rule = "RNA counts > 0",
    scope = "broad_cell_type_condition_union",
    cell_type = label,
    condition_levels = detection$condition_levels,
    n_candidate_genes = length(candidate_genes),
    n_retained_genes = length(retained),
    max_detection = detection$max_detection,
    condition_detection = detection$condition_detection
  )
  message(
    "Stage 1 metabolic target filter | cell_type=", label,
    ";rule=any_condition_detection>=",
    format(threshold, digits = 3, trim = TRUE),
    ";positive=count>0;candidates=", length(candidate_genes),
    ";retained=", length(retained)
  )
  filtered
}

.rc_stage1_filter_pando_detection_features <- function(
    object, pando_motif_args, condition_col,
    rna_assay = "RNA", atac_assay = "ATAC", cell_type = NULL,
    tf_threshold = .rc_stage1_tf_detection_threshold,
    peak_threshold = .rc_stage1_peak_detection_threshold) {
  if (!inherits(object, "Seurat")) {
    stop("Stage 1 regulatory feature filtering requires a Seurat object.",
         call. = FALSE)
  }
  if (!is.list(pando_motif_args)) {
    stop("Stage 1 regulatory feature filtering requires `pando_motif_args`.",
         call. = FALSE)
  }
  motif_tfs <- pando_motif_args$motif_tfs
  if (!is.data.frame(motif_tfs) || ncol(motif_tfs) < 2L) {
    stop(
      "Stage 1 regulatory feature filtering requires materialized ",
      "`pando_motif_args$motif_tfs` with motif and TF columns.",
      call. = FALSE
    )
  }
  tf_threshold <- suppressWarnings(as.numeric(tf_threshold))
  peak_threshold <- suppressWarnings(as.numeric(peak_threshold))
  if (length(tf_threshold) != 1L || !is.finite(tf_threshold) ||
      tf_threshold < 0 || tf_threshold > 1 ||
      length(peak_threshold) != 1L || !is.finite(peak_threshold) ||
      peak_threshold < 0 || peak_threshold > 1) {
    stop("Stage 1 TF/peak detection thresholds must be in [0, 1].",
         call. = FALSE)
  }

  rna_counts <- .rc_get_assay_counts(object, rna_assay)
  atac_counts <- .rc_get_assay_counts(object, atac_assay)
  candidate_tfs <- unique(as.character(motif_tfs[[2L]]))
  candidate_tfs <- candidate_tfs[
    !is.na(candidate_tfs) & nzchar(candidate_tfs) &
      candidate_tfs %in% rownames(rna_counts)
  ]
  if (!length(candidate_tfs)) {
    stop("No motif-linked TF is present in the RNA count matrix.",
         call. = FALSE)
  }
  tf_detection <- .rc_stage1_condition_detection(
    counts = rna_counts, metadata = object@meta.data,
    condition_col = condition_col, features = candidate_tfs,
    keep_matrix = FALSE
  )
  retained_tfs <- candidate_tfs[
    tf_detection$max_detection >= tf_threshold
  ]
  if (!length(retained_tfs)) {
    stop(
      "No motif-linked TF reaches RNA detection >= ",
      format(100 * tf_threshold, trim = TRUE),
      "% (counts > 0) in any retained condition.", call. = FALSE
    )
  }
  filtered_motif_tfs <- motif_tfs[
    as.character(motif_tfs[[2L]]) %in% retained_tfs, , drop = FALSE
  ]

  candidate_peaks <- rownames(atac_counts)
  peak_detection <- .rc_stage1_condition_detection(
    counts = atac_counts, metadata = object@meta.data,
    condition_col = condition_col, features = candidate_peaks,
    keep_matrix = FALSE
  )
  retained_peaks <- candidate_peaks[
    peak_detection$max_detection >= peak_threshold
  ]
  if (!length(retained_peaks)) {
    stop(
      "No ATAC peak reaches detection >= ",
      format(100 * peak_threshold, trim = TRUE),
      "% (counts > 0) in any retained condition.", call. = FALSE
    )
  }
  if (length(retained_peaks) < length(candidate_peaks)) {
    object[[atac_assay]] <- subset(
      object[[atac_assay]], features = retained_peaks
    )
  }

  label <- if (is.null(cell_type) || !length(cell_type) ||
      is.na(cell_type[[1L]]) || !nzchar(trimws(as.character(cell_type[[1L]])))) {
    "<unspecified>"
  } else {
    trimws(as.character(cell_type[[1L]]))
  }
  diagnostics <- list(
    scope = "broad_cell_type_condition_union",
    cell_type = label,
    condition_levels = tf_detection$condition_levels,
    tf_threshold = tf_threshold,
    tf_positive_rule = "RNA counts > 0",
    n_candidate_tfs = length(candidate_tfs),
    n_retained_tfs = length(retained_tfs),
    n_input_motif_tf_rows = nrow(motif_tfs),
    n_retained_motif_tf_rows = nrow(filtered_motif_tfs),
    peak_threshold = peak_threshold,
    peak_positive_rule = "ATAC counts > 0",
    n_candidate_peaks = length(candidate_peaks),
    n_retained_peaks = length(retained_peaks)
  )
  object@misc$regcompass_stage1_regulatory_detection_filter <- diagnostics
  pando_motif_args$motif_tfs <- filtered_motif_tfs
  message(
    "Stage 1 regulatory feature filter | cell_type=", label,
    ";tf_any_condition_detection>=",
    format(tf_threshold, digits = 3, trim = TRUE),
    ";TFs=", length(retained_tfs), "/", length(candidate_tfs),
    ";peak_any_condition_detection>=",
    format(peak_threshold, digits = 3, trim = TRUE),
    ";peaks=", length(retained_peaks), "/", length(candidate_peaks)
  )
  list(
    object = object,
    pando_motif_args = pando_motif_args,
    diagnostics = diagnostics
  )
}

.rc_pando_infer_arg_catalog <- function() {
  list(
    shared = c("tf_cor", "peak_cor", "adjust_method", "padj_threshold"),
    condition = c(
      "rank_action", "min_residual_df", "reference_condition",
      "rna_layer", "peak_layer", "peak_value_type",
      "checkpoint_dir", "resume"
    ),
    standard = c(
      "peak_to_gene_method", "upstream", "downstream", "extend",
      "only_tss", "method", "alpha", "family", "interaction_term",
      "scale", "aggregate_rna_col", "aggregate_peaks_col", "verbose",
      "maxit", "epsilon", "control", "nlambda", "lambda",
      "lambda.min.ratio", "standardize", "nfolds", "type.measure",
      "solver", "bagging_number", "n_jobs", "p_method", "prior",
      "chains", "cores", "iter", "seed", "params", "nrounds",
      "nthread", "ridge_control", "rank_action", "min_residual_df"
    )
  )
}

.rc_route_pando_infer_args <- function(
    args, condition_types = character(), standard_types = character()) {
  if (!is.list(args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  if (is.null(names(args))) names(args) <- rep("", length(args))
  if (any(!nzchar(names(args)))) {
    stop("Every `pando_infer_args` entry must be named.", call. = FALSE)
  }

  if (length(condition_types)) {
    if (!requireNamespace("Pando", quietly = TRUE)) {
      stop("Package 'Pando' is required for condition-GRN Stage 1.",
           call. = FALSE)
    }
    required_api <- c(
      "infer_condition_grn", "condition_grn_fit", "initiate_grn",
      "find_motifs", "LayerData"
    )
    missing_api <- setdiff(required_api, getNamespaceExports("Pando"))
    if (length(missing_api)) {
      stop(
        "Installed Pando lacks required RegCompass condition-GRN API(s): ",
        paste(missing_api, collapse = ", "), ".", call. = FALSE
      )
    }
  }

  obsolete <- intersect(
    names(args),
    c(
      "condition_ridge_control", "condition_e_control", "scheme_e_z", "z",
      "fusion_ratio", "lambda_grid", "lambda_rule", "cv_folds"
    )
  )
  if (length(condition_types) && length(obsolete)) {
    stop(
      "Removed conditional Pando control(s): ",
      paste(obsolete, collapse = ", "),
      ". RegCompass conditional GRNs use fixed production E-star z=0.25 ",
      "coefficients with separate no-fusion exact-edge inference; no ",
      "conditional CV or sensitivity control is exposed.",
      call. = FALSE
    )
  }

  canonical_layers <- list(
    rna_layer = "data", peak_layer = "data", peak_value_type = "normalized"
  )
  supplied_layers <- intersect(names(args), names(canonical_layers))
  if (length(condition_types) && length(supplied_layers)) {
    invalid <- vapply(supplied_layers, function(field) {
      value <- args[[field]]
      !is.character(value) || length(value) != 1L || is.na(value) ||
        !identical(as.character(value), canonical_layers[[field]])
    }, logical(1))
    if (any(invalid)) {
      stop(
        "Condition-GRN Stage 1 requires rna_layer='data', peak_layer='data', ",
        "and peak_value_type='normalized'. Unsupported override(s): ",
        paste(supplied_layers[invalid], collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  catalog <- .rc_pando_infer_arg_catalog()
  condition_allowed <- unique(c(catalog$shared, catalog$condition))
  standard_allowed <- unique(c(catalog$shared, catalog$standard))
  known <- unique(c(condition_allowed, standard_allowed))
  unknown <- setdiff(names(args), known)
  if (length(unknown)) {
    stop(
      "Unsupported `pando_infer_args`: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  condition_args <- args[intersect(names(args), condition_allowed)]
  condition_args <- utils::modifyList(list(
    tf_cor = 0.05, peak_cor = 0.05, adjust_method = "BH",
    padj_threshold = 0.05, rank_action = "mark", min_residual_df = 1L,
    reference_condition = NULL,
    rna_layer = "data", peak_layer = "data", peak_value_type = "normalized",
    checkpoint_dir = NULL, resume = TRUE
  ), condition_args)
  condition_threshold <- suppressWarnings(as.numeric(
    condition_args$padj_threshold
  ))
  if (length(condition_types) &&
      (!identical(toupper(as.character(condition_args$adjust_method)), "BH") ||
       length(condition_threshold) != 1L || !is.finite(condition_threshold) ||
       condition_threshold <= 0 || condition_threshold >= 1)) {
    stop(
      "Canonical RegCompass condition fits require BH adjustment and ",
      "padj_threshold in (0, 1).", call. = FALSE
    )
  }
  condition_args$padj_threshold <- condition_threshold
  if (!is.null(condition_args$checkpoint_dir) &&
      (!is.character(condition_args$checkpoint_dir) ||
       length(condition_args$checkpoint_dir) != 1L ||
       is.na(condition_args$checkpoint_dir) ||
       !nzchar(trimws(condition_args$checkpoint_dir)))) {
    stop("Conditional `checkpoint_dir` must be NULL or one non-empty path.",
         call. = FALSE)
  }
  if (!is.logical(condition_args$resume) ||
      length(condition_args$resume) != 1L || is.na(condition_args$resume)) {
    stop("Conditional `resume` must be either TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(condition_args$reference_condition)) {
    reference <- as.character(condition_args$reference_condition)
    if (length(reference) != 1L || is.na(reference) ||
        !nzchar(trimws(reference)) || reference != trimws(reference)) {
      stop(
        "Conditional `reference_condition` must be NULL or one complete ",
        "predefined condition label.", call. = FALSE
      )
    }
    condition_args$reference_condition <- reference
  }

  standard_args <- args[intersect(names(args), standard_allowed)]
  standard_args <- utils::modifyList(list(
    tf_cor = 0.05, peak_cor = 0.05,
    adjust_method = "BH", padj_threshold = 0.05
  ), standard_args)
  standard_threshold <- suppressWarnings(as.numeric(
    standard_args$padj_threshold
  ))
  if (length(standard_types) &&
      (length(standard_threshold) != 1L || !is.finite(standard_threshold) ||
       standard_threshold <= 0 || standard_threshold >= 1)) {
    stop("Standard Pando `padj_threshold` must be one value in (0, 1).",
         call. = FALSE)
  }
  standard_args$padj_threshold <- standard_threshold
  standard_method <- as.character(standard_args$method %||% "ridge")
  if (length(standard_method) != 1L || is.na(standard_method) ||
      !nzchar(standard_method)) {
    stop("Standard Pando `method` must be one non-empty value.", call. = FALSE)
  }
  standard_args$method <- standard_method
  if (identical(standard_method, "ridge")) {
    standard_args$ridge_control <- standard_args$ridge_control %||% list()
    if (!is.list(standard_args$ridge_control)) {
      stop("Standard Pando `ridge_control` must be a list.", call. = FALSE)
    }
    standard_args$rank_action <- standard_args$rank_action %||% "mark"
    standard_args$min_residual_df <- standard_args$min_residual_df %||% 1L
  } else {
    if ("ridge_control" %in% names(args) && length(args$ridge_control)) {
      stop(
        "`pando_infer_args$ridge_control` requires standard Pando ",
        "`method = \"ridge\"`.", call. = FALSE
      )
    }
    standard_args[c("ridge_control", "rank_action", "min_residual_df")] <- NULL
  }

  condition_only <- setdiff(catalog$condition, catalog$standard)
  standard_only <- setdiff(catalog$standard, catalog$condition)
  disabled_condition <- if (length(standard_types)) {
    intersect(names(args), condition_only)
  } else character()
  disabled_standard <- if (length(condition_types)) {
    intersect(names(args), standard_only)
  } else character()
  diagnostics <- rbind(
    if (length(disabled_condition)) data.frame(
      route = "standard_pando", argument = disabled_condition,
      action = "disabled_condition_only", stringsAsFactors = FALSE
    ),
    if (length(disabled_standard)) data.frame(
      route = "condition_grn", argument = disabled_standard,
      action = "disabled_standard_only", stringsAsFactors = FALSE
    )
  )
  if (is.null(diagnostics)) diagnostics <- data.frame(
    route = character(), argument = character(), action = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(diagnostics)) {
    message(
      "Stage 1 routed `pando_infer_args` by analysis mode; disabled: ",
      paste0(diagnostics$route, "{", diagnostics$argument, "}",
             collapse = ", ")
    )
  }
  list(
    condition = condition_args,
    standard = standard_args,
    diagnostics = diagnostics
  )
}

.rc_standard_pando_single_process_args <- function(args) {
  if (!is.list(args)) return(args)
  controls <- intersect(names(args), c("n_jobs", "cores", "nthread"))
  requested <- if (length(controls)) args[controls] else list()
  for (name in controls) args[[name]] <- 1L
  if (is.list(args$params) && "nthread" %in% names(args$params)) {
    requested$params_nthread <- args$params$nthread
    args$params$nthread <- 1L
  }
  attr(args, "regcompass_pando_parallel_contract") <- list(
    scope = "standard_pando_broad_cell_type_jobs",
    infer_grn_parallel = FALSE,
    inner_worker_limit = 1L,
    overridden_controls = requested
  )
  args
}

.rc_run_standard_pando_celltype_job <- function(
    job, base, extra_args, standard_infer_args,
    outer_parallel, progress_monitor, PANDO_BPPARAM = NULL) {
  if (!is.list(job) || !inherits(job$object, "Seurat")) {
    stop("Invalid standard-Pando cell-type job.", call. = FALSE)
  }
  thread_state <- .rc_set_internal_single_thread()
  on.exit({
    .rc_restore_internal_threads(thread_state)
    invisible(gc(verbose = FALSE, full = TRUE))
  }, add = TRUE)

  job_extra <- extra_args
  motif_args <- job_extra$pando_motif_args %||% list()
  if (outer_parallel && is.list(motif_args) &&
      !is.null(motif_args$cache_dir)) {
    motif_args$cache_dir <- file.path(
      motif_args$cache_dir, .rc_safe_path_component(job$cell_type)
    )
  }
  detection <- .rc_stage1_filter_pando_detection_features(
    object = job$object,
    pando_motif_args = motif_args,
    condition_col = base$condition_col,
    rna_assay = base$rna_assay,
    atac_assay = base$atac_assay,
    cell_type = job$cell_type
  )
  job$object <- detection$object
  job_extra$pando_motif_args <- detection$pando_motif_args
  filtered_gem <- .rc_stage1_filter_gem_metabolic_genes(
    object = job$object, gem = base$gem,
    condition_col = base$condition_col, rna_assay = base$rna_assay,
    cell_type = job$cell_type
  )
  args <- c(base[setdiff(names(base), names(job_extra))], job_extra)
  args$gem <- filtered_gem
  args$object <- job$object
  args$cell_type <- job$cell_type
  args$progress_monitor <- if (outer_parallel) NULL else progress_monitor
  args$outdir <- file.path(
    base$outdir, "standard", .rc_safe_path_component(job$cell_type)
  )
  fixed_infer_args <- .rc_standard_pando_single_process_args(
    standard_infer_args
  )
  pando_parallel_contract <- attr(
    fixed_infer_args, "regcompass_pando_parallel_contract", exact = TRUE
  )
  attr(fixed_infer_args, "regcompass_pando_parallel_contract") <- NULL
  method <- as.character(fixed_infer_args$method %||% "ridge")
  ridge_inner_parallel <- identical(method, "ridge") && !isTRUE(outer_parallel) &&
    !is.null(PANDO_BPPARAM) && !identical(PANDO_BPPARAM, FALSE) &&
    .rc_bpparam_worker_limit(PANDO_BPPARAM, default = 1L) > 1L
  args$pando_infer_args <- fixed_infer_args
  args$parallel <- ridge_inner_parallel
  args$BPPARAM <- if (ridge_inner_parallel) PANDO_BPPARAM else NULL
  if (is.list(pando_parallel_contract)) {
    pando_parallel_contract$scope <- if (ridge_inner_parallel) {
      "standard_pando_single_celltype_target_parallel"
    } else {
      "standard_pando_celltype_parallel_or_serial"
    }
    pando_parallel_contract$infer_grn_parallel <- ridge_inner_parallel
    pando_parallel_contract$inner_worker_limit <- if (ridge_inner_parallel) {
      .rc_bpparam_worker_limit(PANDO_BPPARAM, default = 1L)
    } else 1L
    pando_parallel_contract$nested_parallel <- FALSE
  }
  value <- do.call(.rc_fit_standard_pando_by_cell_type, args)
  if (is.list(value$normalization_policy)) {
    value$normalization_policy$parallel_contract <- pando_parallel_contract
    value$normalization_policy$metabolic_target_filter <-
      "RNA counts > 0 in >=20% of cells in any retained condition within the broad cell type"
    value$normalization_policy$regulatory_feature_detection_filter <-
      detection$diagnostics
  }
  list(cell_type = job$cell_type, route = "standard_pando", result = value)
}

.rc_run_condition_pando_batch <- function(
    object, condition_types, base, extra_args, condition_infer_args,
    parallel, BPPARAM, progress_monitor) {
  if (!length(condition_types)) return(NULL)
  worker_limit <- if (isTRUE(parallel) && !is.null(BPPARAM) &&
      !identical(BPPARAM, FALSE)) {
    .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  } else 1L
  values <- vector("list", length(condition_types))
  names(values) <- condition_types
  for (i in seq_along(condition_types)) {
    type <- condition_types[[i]]
    cells <- rownames(object@meta.data)[
      as.character(object@meta.data[[base$celltype_col]]) == type
    ]
    if (!length(cells)) {
      stop("No cells remain for condition-GRN cell type `", type, "`.",
           call. = FALSE)
    }
    one <- subset(object, cells = cells)
    one <- .rc_stage1_pando_working_object(
      one, rna_assay = base$rna_assay, atac_assay = base$atac_assay
    )
    job_extra <- extra_args
    detection <- .rc_stage1_filter_pando_detection_features(
      object = one,
      pando_motif_args = job_extra$pando_motif_args %||% list(),
      condition_col = base$condition_col,
      rna_assay = base$rna_assay,
      atac_assay = base$atac_assay,
      cell_type = type
    )
    one <- detection$object
    job_extra$pando_motif_args <- detection$pando_motif_args
    filtered_gem <- .rc_stage1_filter_gem_metabolic_genes(
      object = one, gem = base$gem,
      condition_col = base$condition_col, rna_assay = base$rna_assay,
      cell_type = type
    )
    args <- c(base[setdiff(names(base), names(job_extra))], job_extra)
    args$gem <- filtered_gem
    args$object <- one
    args$cell_type <- type
    args$outdir <- file.path(
      base$outdir, "condition", .rc_safe_path_component(type)
    )
    args$pando_infer_args <- condition_infer_args
    args$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE
    args$progress_monitor <- progress_monitor
    values[[i]] <- do.call(.rc_fit_condition_grns_by_cell_type, args)
    if (is.list(values[[i]]$normalization_policy)) {
      values[[i]]$normalization_policy$metabolic_target_filter <-
        "RNA counts > 0 in >=20% of cells in any retained condition within the broad cell type"
      values[[i]]$normalization_policy$regulatory_feature_detection_filter <-
        detection$diagnostics
    }
    one <- NULL
    detection <- NULL
    filtered_gem <- NULL
    job_extra <- NULL
    args <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
  }
  answer <- if (length(values) == 1L) {
    values[[1L]]
  } else {
    .rc_merge_condition_job_results(values)
  }
  plan <- list(
    scope = "cell_type_serial_target_parallel",
    cell_types = condition_types,
    workers = worker_limit,
    outer_celltype_parallel = FALSE,
    inner_target_parallel = isTRUE(parallel) && worker_limit > 1L,
    nested_parallel = FALSE,
    worker_budget_shared_sequentially = TRUE,
    memory_policy = paste(
      "one cell type resident at a time; Pando target workers receive",
      "target-specific compact payloads"
    )
  )
  if (!is.list(answer$pando_execution_summary)) {
    answer$pando_execution_summary <- list()
  }
  answer$pando_execution_summary$parallel_plan <- plan
  if (!is.list(answer$normalization_policy)) {
    answer$normalization_policy <- list()
  }
  answer$normalization_policy$parallel_contract <- plan
  answer$normalization_policy$metabolic_target_filter <-
    "RNA counts > 0 in >=20% of cells in any retained condition within the broad cell type"
  answer$normalization_policy$regulatory_feature_detection_filter <- list(
    scope = "broad_cell_type_condition_union",
    tf = "RNA counts > 0 in >=5% of cells in any retained condition",
    peak = "ATAC counts > 0 in >=5% of cells in any retained condition"
  )
  answer
}

.rc_fit_pando_by_celltype_route <- function(
    object, gem, outdir, genome, pfm, species, condition_col, celltype_col,
    condition_types, standard_types, rna_assay, atac_assay,
    extra_args, condition_infer_args, standard_infer_args,
    parallel, BPPARAM, progress_monitor) {
  if (!length(condition_types) && !length(standard_types)) {
    stop("No Pando cell-type job was selected.", call. = FALSE)
  }
  base <- list(
    gem = gem, outdir = outdir, genome = genome,
    pfm = pfm, species = species, condition_col = condition_col,
    celltype_col = celltype_col, rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  standard_method <- as.character(standard_infer_args$method %||% "ridge")
  standard_ridge <- identical(standard_method, "ridge")
  worker_limit <- if (isTRUE(parallel) && !is.null(BPPARAM) &&
      !identical(BPPARAM, FALSE)) {
    .rc_bpparam_worker_limit(BPPARAM, default = 1L)
  } else 1L

  .rc_step_monitor_event(
    progress_monitor, "cell_type_execution_plan",
    paste(
      "condition GRNs use pooled/global plus condition exact union,",
      "E-star z=0.25 production plus separated exact-edge inference;",
      "standard Pando uses", standard_method,
      if (standard_ridge) "through the independent K=1 ridge route" else ""
    ),
    current = 5L,
    context = list(
      condition_cell_types = length(condition_types),
      standard_cell_types = length(standard_types),
      condition_reference = condition_infer_args$reference_condition %||%
        "<first-retained>",
      condition_parallel_scope = if (length(condition_types) &&
          isTRUE(parallel) && worker_limit > 1L) "target" else "serial",
      standard_parallel_scope = if (length(standard_types) && standard_ridge &&
          isTRUE(parallel) && worker_limit > 1L) {
        "target"
      } else if (length(standard_types) > 1L && !standard_ridge &&
                 isTRUE(parallel)) {
        "cell_type"
      } else {
        "serial"
      },
      nested_parallel = FALSE,
      memory_policy = paste(
        "single parallel level; one ridge/E-star cell type resident at a time;",
        "target-specific Pando worker payloads"
      )
    )
  )

  condition_result <- .rc_run_condition_pando_batch(
    object = object, condition_types = condition_types, base = base,
    extra_args = extra_args, condition_infer_args = condition_infer_args,
    parallel = parallel, BPPARAM = BPPARAM,
    progress_monitor = progress_monitor
  )
  invisible(gc(verbose = FALSE, full = TRUE))

  standard_values <- list()
  standard_outer_parallel <- FALSE
  if (length(standard_types)) {
    if (standard_ridge) {
      for (i in seq_along(standard_types)) {
        type <- standard_types[[i]]
        cells <- rownames(object@meta.data)[
          as.character(object@meta.data[[celltype_col]]) == type
        ]
        one <- subset(object, cells = cells)
        one <- .rc_stage1_pando_working_object(
          one, rna_assay = rna_assay, atac_assay = atac_assay
        )
        job <- list(cell_type = type, object = one)
        executed <- .rc_run_standard_pando_celltype_job(
          job = job, base = base, extra_args = extra_args,
          standard_infer_args = standard_infer_args,
          outer_parallel = FALSE, progress_monitor = progress_monitor,
          PANDO_BPPARAM = if (isTRUE(parallel)) BPPARAM else NULL
        )
        standard_values[[type]] <- executed$result
        one <- NULL
        job <- NULL
        executed <- NULL
        invisible(gc(verbose = FALSE, full = TRUE))
      }
    } else {
      standard_source <- .rc_stage1_pando_working_object(
        object, rna_assay = rna_assay, atac_assay = atac_assay
      )
      standard_inputs <- lapply(standard_types, function(type) {
        cells <- rownames(standard_source@meta.data)[
          as.character(standard_source@meta.data[[celltype_col]]) == type
        ]
        list(
          cell_type = type,
          object = subset(standard_source, cells = cells)
        )
      })
      standard_outer_parallel <- isTRUE(parallel) &&
        length(standard_inputs) > 1L
      executed <- rc_parallel_lapply(
        standard_inputs, .rc_run_standard_pando_celltype_job,
        BPPARAM = if (standard_outer_parallel) BPPARAM else FALSE,
        base = base, extra_args = extra_args,
        standard_infer_args = standard_infer_args,
        outer_parallel = standard_outer_parallel,
        progress_monitor = progress_monitor, PANDO_BPPARAM = NULL
      )
      standard_values <- lapply(executed, `[[`, "result")
      names(standard_values) <- vapply(
        executed, `[[`, character(1), "cell_type"
      )
      standard_source <- NULL
      standard_inputs <- NULL
      executed <- NULL
      invisible(gc(verbose = FALSE, full = TRUE))
    }
  }

  answer <- .rc_merge_pando_results(
    condition_result = condition_result,
    standard_results = standard_values,
    condition_types = condition_types, standard_types = standard_types,
    condition_col = condition_col, celltype_col = celltype_col,
    outdir = outdir
  )
  condition_plan <- if (!is.null(condition_result) &&
      is.list(condition_result$pando_execution_summary)) {
    condition_result$pando_execution_summary$parallel_plan %||% list()
  } else {
    list()
  }
  answer$pando_execution_plan <- list(
    scope = if (length(condition_types) && length(standard_types)) {
      "condition_then_standard_target_parallel"
    } else if (length(condition_types)) {
      "condition_target_parallel"
    } else {
      "standard_target_parallel"
    },
    condition_cell_types = condition_types,
    standard_cell_types = standard_types,
    condition_parallel_plan = condition_plan,
    standard_method = standard_method,
    standard_outer_parallel = standard_outer_parallel,
    standard_target_parallel = standard_ridge && isTRUE(parallel) &&
      worker_limit > 1L,
    nested_parallel = FALSE,
    worker_budget_shared_sequentially = TRUE,
    memory_policy = paste(
      "one ridge/E-star cell type resident at a time; target-specific Pando",
      "payloads; worker and batch temporaries released after completion"
    )
  )
  if (!is.list(answer$normalization_policy)) {
    answer$normalization_policy <- list()
  }
  answer$normalization_policy$stage1_detection_filters <- list(
    scope = "broad_cell_type_condition_union",
    positive = "raw assay counts > 0",
    metabolic_target_rna = ">=20% in any retained condition",
    tf_rna = ">=5% in any retained condition",
    peak_atac = ">=5% in any retained condition",
    application = paste(
      "metabolic target genes restrict the Pando target universe; TF filtering",
      "restricts motif-to-TF mappings; peak filtering restricts the ATAC feature",
      "space before Pando candidate discovery"
    )
  )
  answer
}
