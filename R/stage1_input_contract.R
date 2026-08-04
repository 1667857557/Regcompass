.rc_stage1_min_cells_fixed <- 300L

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
  is.function(method) &&
    all(c("cache_dir", "reuse_cache") %in% names(formals(method)))
}

.rc_resolve_stage1_min_cells_contract <- function(pando_args) {
  if (!is.list(pando_args)) {
    stop("`pando_args` must be a list.", call. = FALSE)
  }
  supplied <- pando_args$min_cells %||% .rc_stage1_min_cells_fixed
  supplied <- suppressWarnings(as.integer(supplied[[1L]]))
  if (!is.finite(supplied) || supplied != .rc_stage1_min_cells_fixed) {
    message("Stage 1 `min_cells` is fixed at 300; overriding the supplied value.")
  }
  pando_args$min_cells <- .rc_stage1_min_cells_fixed
  list(min_cells = .rc_stage1_min_cells_fixed, pando_args = pando_args)
}

.rc_filter_stage1_groups_by_min_cells <- function(
    object, condition_col, celltype_col, cell_type, min_cells) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  .rc_validate_celltype_metadata(object@meta.data, celltype_col)
  observed_type <- trimws(as.character(object@meta.data[[celltype_col]]))
  available_type <- unique(observed_type)
  requested_type <- if (is.null(cell_type)) {
    available_type
  } else {
    unique(trimws(as.character(cell_type)))
  }
  missing_type <- setdiff(requested_type, available_type)
  if (length(missing_type)) {
    stop("Requested cell types were not found: ",
         paste(missing_type, collapse = ", "), call. = FALSE)
  }
  selected <- observed_type %in% requested_type
  if (!any(selected)) {
    stop("No cells remain after applying `cell_type`.", call. = FALSE)
  }

  has_condition <- is.character(condition_col) && length(condition_col) == 1L &&
    !is.na(condition_col) && nzchar(condition_col) &&
    condition_col %in% colnames(object@meta.data)
  if (has_condition) {
    observed_condition <- trimws(as.character(object@meta.data[[condition_col]]))
    invalid_condition <- selected &
      (is.na(object@meta.data[[condition_col]]) | !nzchar(observed_condition))
    if (any(invalid_condition)) {
      stop("Condition metadata contain missing or empty values in requested cell types.",
           call. = FALSE)
    }
    condition_levels <- unique(observed_condition[selected])
  } else {
    observed_condition <- rep(NA_character_, nrow(object@meta.data))
    condition_levels <- character()
  }

  condition_design <- length(condition_levels) >= 2L
  if (condition_design) {
    count_matrix <- table(
      factor(observed_type[selected], levels = requested_type),
      factor(observed_condition[selected], levels = condition_levels)
    )
    retained_stratum <- count_matrix >= as.integer(min_cells)
    retained_condition_count <- base::rowSums(retained_stratum)
    retained_type <- rownames(count_matrix)[retained_condition_count >= 1L]
    condition_type <- rownames(count_matrix)[retained_condition_count >= 2L]
    standard_type <- rownames(count_matrix)[retained_condition_count == 1L]

    diagnostics <- do.call(rbind, lapply(requested_type, function(type) {
      stratum_ok <- as.logical(retained_stratum[type, condition_levels])
      type_mode <- if (type %in% condition_type) {
        "condition_grn"
      } else if (type %in% standard_type) {
        "standard_pando"
      } else {
        "excluded"
      }
      data.frame(
        cell_type = type,
        condition = condition_levels,
        n_cells = as.integer(count_matrix[type, condition_levels]),
        retained_stratum = stratum_ok,
        retained_cell_type = type %in% retained_type,
        cell_type_analysis_mode = type_mode,
        retained = stratum_ok,
        fit_status = ifelse(
          !stratum_ok,
          "excluded_below_min_cells",
          ifelse(type %in% condition_type,
                 "eligible_condition_pando",
                 "eligible_standard_pando_single_condition")
        ),
        n_retained_conditions = as.integer(retained_condition_count[[type]]),
        threshold = as.integer(min_cells),
        threshold_scope = "condition_x_cell_type",
        stringsAsFactors = FALSE
      )
    }))

    dropped <- diagnostics[!diagnostics$retained_stratum, , drop = FALSE]
    if (nrow(dropped)) {
      message(
        "Stage 1 excluded condition x cell-type strata below min_cells=300: ",
        paste0(dropped$cell_type, "{", dropped$condition, "}=",
               dropped$n_cells, collapse = "; ")
      )
    }
    if (!length(retained_type)) {
      stop("No condition x cell-type stratum reaches min_cells=300.",
           call. = FALSE)
    }
    if (length(standard_type)) {
      message(
        "Using standard Pando for cell types with one retained condition: ",
        paste(standard_type, collapse = ", ")
      )
    }

    keep_cells <- rep(FALSE, nrow(object@meta.data))
    selected_rows <- which(selected)
    type_index <- match(observed_type[selected_rows], rownames(retained_stratum))
    condition_index <- match(
      observed_condition[selected_rows], colnames(retained_stratum)
    )
    valid <- !is.na(type_index) & !is.na(condition_index)
    selected_keep <- rep(FALSE, length(selected_rows))
    selected_keep[valid] <- retained_stratum[cbind(
      type_index[valid], condition_index[valid]
    )]
    keep_cells[selected_rows] <- selected_keep
    analysis_mode <- if (length(condition_type) && length(standard_type)) {
      "mixed_pando"
    } else if (length(condition_type)) {
      "condition_grn"
    } else {
      "standard_pando"
    }
  } else {
    counts <- table(factor(observed_type[selected], levels = requested_type))
    retained_type <- names(counts)[as.integer(counts) >= min_cells]
    condition_type <- character()
    standard_type <- retained_type
    diagnostics <- data.frame(
      cell_type = names(counts),
      condition = if (length(condition_levels) == 1L) {
        condition_levels[[1L]]
      } else {
        NA_character_
      },
      n_cells = as.integer(counts),
      retained_stratum = as.integer(counts) >= min_cells,
      retained_cell_type = names(counts) %in% retained_type,
      cell_type_analysis_mode = ifelse(
        names(counts) %in% retained_type, "standard_pando", "excluded"
      ),
      retained = names(counts) %in% retained_type,
      fit_status = ifelse(
        names(counts) %in% retained_type,
        "eligible_standard_pando",
        "excluded_below_min_cells"
      ),
      n_retained_conditions = if (length(condition_levels) == 1L) {
        as.integer(names(counts) %in% retained_type)
      } else {
        NA_integer_
      },
      threshold = as.integer(min_cells),
      threshold_scope = "cell_type",
      stringsAsFactors = FALSE
    )
    if (!length(retained_type)) {
      stop("No requested cell type reaches min_cells=300.", call. = FALSE)
    }
    keep_cells <- selected & observed_type %in% retained_type
    analysis_mode <- "standard_pando"
  }

  retained_cells <- rownames(object@meta.data)[keep_cells]
  filtered <- subset(object, cells = retained_cells)
  filtered@misc$regcompass_stage1_group_filter <- diagnostics
  list(
    object = filtered,
    retained_cell_types = retained_type,
    condition_pando_cell_types = condition_type,
    standard_pando_cell_types = standard_type,
    diagnostics = diagnostics,
    analysis_mode = analysis_mode,
    condition_levels = if (condition_design) {
      unique(observed_condition[keep_cells])
    } else {
      condition_levels
    }
  )
}

.rc_build_stage_analysis_cell_set <- function(
    object, condition_col, celltype_col, cell_type, pando_args) {
  threshold <- .rc_resolve_stage1_min_cells_contract(pando_args)
  groups <- .rc_filter_stage1_groups_by_min_cells(
    object = object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    min_cells = threshold$min_cells
  )
  analysis_types <- unique(as.character(groups$retained_cell_types))
  observed_type <- trimws(as.character(groups$object@meta.data[[celltype_col]]))
  analysis_cells <- rownames(groups$object@meta.data)[
    observed_type %in% analysis_types
  ]
  if (!length(analysis_cells)) {
    stop("No cells remain for Stage 1 Pando analysis.", call. = FALSE)
  }
  analysis_object <- subset(groups$object, cells = analysis_cells)
  analysis_cells <- unname(as.character(colnames(analysis_object)))
  list(
    object = analysis_object,
    retained_cells = analysis_cells,
    retained_cell_types = analysis_types,
    condition_pando_cell_types =
      unique(as.character(groups$condition_pando_cell_types)),
    standard_pando_cell_types =
      unique(as.character(groups$standard_pando_cell_types)),
    diagnostics = groups$diagnostics,
    analysis_mode = groups$analysis_mode,
    condition_levels = groups$condition_levels,
    min_cells = threshold$min_cells,
    pando_args = threshold$pando_args
  )
}

.rc_validate_stage1_cell_set <- function(grn) {
  if (!inherits(grn, "regcompass_grn_step")) {
    stop("`grn` must be the output of `rc_regcompass_step_grn()`.",
         call. = FALSE)
  }
  contract <- grn$cell_filter
  if (!is.list(contract)) {
    stop(
      "The Stage 1 result has no cell-set contract. Rerun Stage 1 with the current RegCompassR version.",
      call. = FALSE
    )
  }
  cells <- unname(as.character(contract$retained_cells))
  types <- unname(as.character(contract$retained_cell_types))
  if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) ||
      anyDuplicated(cells)) {
    stop("The Stage 1 retained-cell contract is invalid.", call. = FALSE)
  }
  if (!length(types) || anyNA(types) || any(!nzchar(types))) {
    stop("The Stage 1 retained-cell-type contract is invalid.", call. = FALSE)
  }
  contract$retained_cells <- cells
  contract$retained_cell_types <- unique(types)
  contract
}

.rc_subset_to_stage1_cell_set <- function(object, contract) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  expected <- unname(as.character(contract$retained_cells))
  available <- unname(as.character(colnames(object)))
  missing <- setdiff(expected, available)
  if (length(missing)) {
    stop(
      "Stage 2 input is missing ", length(missing),
      " cell(s) retained by Stage 1; first missing ID: ", missing[[1L]],
      call. = FALSE
    )
  }
  filtered <- subset(object, cells = expected)
  observed <- unname(as.character(colnames(filtered)))
  exact_set <- length(observed) == length(expected) &&
    !anyDuplicated(observed) && setequal(observed, expected)
  if (!exact_set) {
    stop("Stage 2 could not reproduce the exact Stage 1 cell set.",
         call. = FALSE)
  }
  order_matches_stage1 <- identical(observed, expected)
  filtered@misc$regcompass_cross_stage_cell_set <- list(
    source = contract$source %||% "stage1_grn_result",
    n_cells = length(expected),
    retained_cell_types = contract$retained_cell_types,
    min_cells = contract$min_cells %||% .rc_stage1_min_cells_fixed,
    exact_cell_set = TRUE,
    order_matches_stage1 = order_matches_stage1,
    order_policy = if (order_matches_stage1) {
      "stage1_contract_order"
    } else {
      "seurat_native_order_exact_set"
    },
    stage1_order_index = match(observed, expected)
  )
  filtered
}
