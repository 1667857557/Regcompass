# Final Stage 1 threshold contract. This file is collated after the workflow
# hardening overrides so one fixed `min_cells` value drives prefiltering and the
# downstream standard/condition Pando calls.

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
  list(min_cells = pando_args$min_cells, pando_args = pando_args)
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
    stop(
      "Requested cell types were not found: ",
      paste(missing_type, collapse = ", "),
      call. = FALSE
    )
  }
  selected <- observed_type %in% requested_type
  if (!any(selected)) {
    stop("No cells remain after applying `cell_type`.", call. = FALSE)
  }

  has_condition <- is.character(condition_col) && length(condition_col) == 1L &&
    !is.na(condition_col) && nzchar(condition_col) &&
    condition_col %in% colnames(object@meta.data)
  if (has_condition) {
    observed_condition <- trimws(as.character(
      object@meta.data[[condition_col]]
    ))
    invalid_condition <- selected &
      (is.na(object@meta.data[[condition_col]]) | !nzchar(observed_condition))
    if (any(invalid_condition)) {
      stop(
        "Condition metadata contain missing or empty values in requested cell types.",
        call. = FALSE
      )
    }
    condition_levels <- unique(observed_condition[selected])
  } else {
    observed_condition <- rep(NA_character_, nrow(object@meta.data))
    condition_levels <- character()
  }

  condition_mode <- length(condition_levels) >= 2L
  if (condition_mode) {
    count_matrix <- table(
      factor(observed_type[selected], levels = requested_type),
      factor(observed_condition[selected], levels = condition_levels)
    )
    retained_stratum <- count_matrix >= as.integer(min_cells)
    retained_condition_count <- Matrix::rowSums(retained_stratum)
    retained_type <- rownames(count_matrix)[retained_condition_count >= 2L]

    diagnostics <- do.call(rbind, lapply(requested_type, function(type) {
      data.frame(
        cell_type = type,
        condition = condition_levels,
        n_cells = as.integer(count_matrix[type, condition_levels]),
        retained_stratum = as.logical(
          retained_stratum[type, condition_levels]
        ),
        retained_cell_type = type %in% retained_type,
        retained = type %in% retained_type & as.logical(
          retained_stratum[type, condition_levels]
        ),
        n_retained_conditions = as.integer(
          retained_condition_count[[type]]
        ),
        threshold = as.integer(min_cells),
        threshold_scope = "condition_x_cell_type",
        analysis_mode = "condition_grn",
        stringsAsFactors = FALSE
      )
    }))

    dropped_strata <- diagnostics[
      !diagnostics$retained_stratum, , drop = FALSE
    ]
    if (nrow(dropped_strata)) {
      message(
        "Stage 1 excluded condition x cell-type strata below min_cells=300 ",
        "before normalization: ",
        paste0(
          dropped_strata$cell_type, "{", dropped_strata$condition, "}=",
          dropped_strata$n_cells,
          collapse = "; "
        )
      )
    }

    dropped_type <- setdiff(requested_type, retained_type)
    if (length(dropped_type)) {
      message(
        "Stage 1 excluded cell types with fewer than two retained conditions ",
        "after condition-level min_cells filtering: ",
        paste0(
          dropped_type, "=", as.integer(retained_condition_count[dropped_type]),
          " retained condition(s)",
          collapse = "; "
        )
      )
    }

    if (!length(retained_type)) {
      detail <- vapply(requested_type, function(type) {
        paste0(
          type, "{",
          paste0(
            condition_levels, "=", as.integer(count_matrix[type, condition_levels]),
            collapse = ","
          ),
          "}"
        )
      }, character(1))
      stop(
        "No requested cell type retains at least two conditions with ",
        "min_cells=300. Observed: ", paste(detail, collapse = "; "),
        call. = FALSE
      )
    }

    keep_cells <- vapply(seq_len(nrow(object@meta.data)), function(i) {
      type <- observed_type[[i]]
      condition <- observed_condition[[i]]
      selected[[i]] && type %in% retained_type &&
        condition %in% condition_levels &&
        isTRUE(retained_stratum[type, condition])
    }, logical(1))
  } else {
    counts <- table(factor(observed_type[selected], levels = requested_type))
    retained_type <- names(counts)[as.integer(counts) >= min_cells]
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
      retained = names(counts) %in% retained_type,
      n_retained_conditions = if (length(condition_levels) == 1L) {
        as.integer(names(counts) %in% retained_type)
      } else {
        NA_integer_
      },
      threshold = as.integer(min_cells),
      threshold_scope = "cell_type",
      analysis_mode = "standard_pando",
      stringsAsFactors = FALSE
    )
    dropped_type <- setdiff(requested_type, retained_type)
    if (length(dropped_type)) {
      message(
        "Stage 1 excluded cell types with fewer than min_cells=300 before ",
        "normalization: ",
        paste0(dropped_type, "=", as.integer(counts[dropped_type]), collapse = ", ")
      )
    }
    if (!length(retained_type)) {
      stop(
        "No requested cell type reaches min_cells=300.",
        call. = FALSE
      )
    }
    keep_cells <- selected & observed_type %in% retained_type
  }

  retained_cells <- rownames(object@meta.data)[keep_cells]
  filtered <- subset(object, cells = retained_cells)
  filtered@misc$regcompass_stage1_group_filter <- diagnostics
  retained_condition_levels <- if (condition_mode) {
    unique(observed_condition[keep_cells])
  } else {
    condition_levels
  }
  list(
    object = filtered,
    retained_cell_types = retained_type,
    diagnostics = diagnostics,
    analysis_mode = if (condition_mode) "condition_grn" else "standard_pando",
    condition_levels = retained_condition_levels
  )
}

#' Infer regulatory evidence using the fixed Stage 1 `min_cells` contract
#'
#' Stage 1 fixes `pando_args$min_cells` at 300. In standard-Pando mode, each
#' broad cell type must contain at least 300 cells. In condition-Pando mode,
#' every condition-by-cell-type stratum is checked independently: strata below
#' 300 cells are removed while qualifying conditions of the same cell type are
#' retained. A cell type enters condition Pando only when at least two qualifying
#' conditions remain. The same value is then passed unchanged into Pando.
#' Globally zero ATAC peaks are removed before TF-IDF or motif analysis. By
#' default stale Signac fragment references are cleared because Stage 1 uses the
#' in-memory peak matrix and genome sequence rather than fragment files.
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
  threshold_contract <- .rc_resolve_stage1_min_cells_contract(pando_args)
  pando_args <- threshold_contract$pando_args
  min_cells <- threshold_contract$min_cells

  preserve_fragments <- .rc_validate_stage1_fragment_policy(fragment_files)
  filtered_groups <- .rc_filter_stage1_groups_by_min_cells(
    object = object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    min_cells = min_cells
  )
  object <- filtered_groups$object
  filtered_groups$diagnostics$threshold_source <-
    "fixed_pando_args_min_cells"
  object@misc$regcompass_stage1_group_filter <- filtered_groups$diagnostics
  object@misc$regcompass_stage1_min_cells_contract <- list(
    min_cells = min_cells,
    source = "pando_args$min_cells",
    fixed = TRUE,
    analysis_mode = filtered_groups$analysis_mode,
    threshold_scope = if (identical(
      filtered_groups$analysis_mode, "condition_grn"
    )) {
      "condition_x_cell_type_independent_then_minimum_two_conditions"
    } else {
      "cell_type"
    },
    condition_levels = filtered_groups$condition_levels,
    applied_before_normalization = TRUE,
    passed_to_pando = TRUE
  )

  if (!preserve_fragments) {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  }
  zero_filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Stage 1 min_cells prefilter"
  )
  object <- zero_filtered$object
  object@misc$regcompass_stage1_zero_peak_filter <- zero_filtered$diagnostics
  object@misc$regcompass_stage1_fragment_policy <- list(
    fragment_files = preserve_fragments,
    policy = if (preserve_fragments) "preserve" else "clear_before_stage1"
  )

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
    cell_type = if (is.null(cell_type)) NULL else
      filtered_groups$retained_cell_types,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    pando_args = pando_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
}
