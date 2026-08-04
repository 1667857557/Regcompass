# Backward compatibility and downstream validation for Pando condition fits.

.rc_complete_pando_condition_fit_contract <- function(fit) {
  if (!inherits(fit, "ConditionGRNFit")) return(fit)
  field <- "dictionary_preprocessing_provenance_verified"
  if (!field %in% names(fit)) {
    fit[[field]] <- isTRUE(attr(
      fit$edge_dictionary,
      "preprocessing_provenance_verified",
      exact = TRUE
    ))
  }
  fit
}

.rc_complete_pando_condition_fit_contracts <- function(fits) {
  if (inherits(fits, "ConditionGRNFit")) {
    return(.rc_complete_pando_condition_fit_contract(fits))
  }
  if (!is.list(fits)) return(fits)
  lapply(fits, .rc_complete_pando_condition_fit_contract)
}

.rc_require_pando_condition_grn_fit_base <-
  .rc_require_pando_condition_grn_fit

.rc_require_pando_condition_grn_fit <- function(fit) {
  fit <- .rc_complete_pando_condition_fit_contract(fit)
  .rc_require_pando_condition_grn_fit_base(fit)

  downstream_required <- c(
    "condition_col", "cell_type_col", "fit_engine",
    "coefficient_scale", "target_genes"
  )
  if (!all(downstream_required %in% names(fit))) {
    stop(
      "Pando condition fit lacks downstream RegCompass fields: ",
      paste(setdiff(downstream_required, names(fit)), collapse = ", "),
      ".", call. = FALSE
    )
  }
  scalar_text <- function(value) {
    is.character(value) && length(value) == 1L &&
      !is.na(value) && nzchar(trimws(value))
  }
  if (!scalar_text(fit$cell_type) ||
      !scalar_text(fit$condition_col) ||
      !scalar_text(fit$cell_type_col) ||
      identical(fit$condition_col, fit$cell_type_col) ||
      !scalar_text(fit$fit_engine) ||
      !scalar_text(fit$coefficient_scale)) {
    stop("Pando condition fit identifiers and model labels are invalid.",
         call. = FALSE)
  }

  levels <- as.character(fit$condition_levels)
  if (length(levels) < 2L || anyNA(levels) || any(!nzchar(levels)) ||
      anyDuplicated(levels)) {
    stop("Pando condition levels must contain at least two unique labels.",
         call. = FALSE)
  }
  cells_by_condition <- fit$condition_cell_ids
  cell_list_names <- names(cells_by_condition)
  if (!is.list(cells_by_condition) || is.null(cell_list_names) ||
      anyNA(cell_list_names) || any(!nzchar(cell_list_names)) ||
      anyDuplicated(cell_list_names) ||
      !all(levels %in% cell_list_names)) {
    stop("Pando condition cell IDs are not uniquely named for every condition.",
         call. = FALSE)
  }
  cells_by_condition <- cells_by_condition[levels]
  if (any(lengths(cells_by_condition) < 1L)) {
    stop("Every Pando fitted condition must contain at least one cell.",
         call. = FALSE)
  }
  cells <- as.character(unlist(cells_by_condition, use.names = FALSE))
  if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) ||
      anyDuplicated(cells)) {
    stop("Pando fitted cells must be complete and condition-disjoint.",
         call. = FALSE)
  }

  targets <- unique(as.character(fit$target_genes))
  target_key <- toupper(targets)
  if (!length(targets) || anyNA(targets) || any(!nzchar(targets)) ||
      anyDuplicated(target_key)) {
    stop("Pando fitted target genes are empty or case-ambiguous.",
         call. = FALSE)
  }

  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  required_fit <- c("target", "condition", "rsq", "fit_status")
  if (!all(required_fit %in% colnames(fit_table)) || !nrow(fit_table)) {
    stop("Pando target-level fit diagnostics are incomplete.",
         call. = FALSE)
  }
  fit_key <- paste(
    toupper(as.character(fit_table$target)),
    as.character(fit_table$condition), sep = "\001"
  )
  if (anyNA(fit_key) || anyDuplicated(fit_key) ||
      any(!as.character(fit_table$condition) %in% levels)) {
    stop("Pando target-level fit diagnostics are duplicated or mislabelled.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_validate_pando_fit_metadata_frame <- function(
    metadata, fits, condition_col, celltype_col) {
  if (!is.data.frame(metadata) ||
      !all(c(condition_col, celltype_col) %in% colnames(metadata)) ||
      is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
    stop("Pando object metadata cannot validate condition fit cell mappings.",
         call. = FALSE)
  }
  if (inherits(fits, "ConditionGRNFit")) fits <- list(fits)
  if (!is.list(fits) || !length(fits)) {
    stop("No Pando condition fits are available for metadata validation.",
         call. = FALSE)
  }
  for (fit in fits) {
    if (!identical(as.character(fit$condition_col), condition_col) ||
        !identical(as.character(fit$cell_type_col), celltype_col)) {
      stop(
        "Pando fit metadata columns do not match the RegCompass request: ",
        "fit condition_col='", as.character(fit$condition_col),
        "', cell_type_col='", as.character(fit$cell_type_col),
        "'; requested condition_col='", condition_col,
        "', cell_type_col='", celltype_col, "'.",
        call. = FALSE
      )
    }
    levels <- as.character(fit$condition_levels)
    cells_by_condition <- fit$condition_cell_ids[levels]
    for (condition in levels) {
      cells <- as.character(cells_by_condition[[condition]])
      missing <- setdiff(cells, rownames(metadata))
      if (length(missing)) {
        stop(
          "Pando fit references cells absent from its stored object; first ",
          "missing ID: ", missing[[1L]], ".", call. = FALSE
        )
      }
      observed_condition <- as.character(metadata[cells, condition_col])
      observed_celltype <- as.character(metadata[cells, celltype_col])
      if (anyNA(observed_condition) || anyNA(observed_celltype) ||
          any(observed_condition != condition) ||
          any(observed_celltype != as.character(fit$cell_type))) {
        stop(
          "Pando fit cell assignments disagree with stored object metadata ",
          "for cell type '", as.character(fit$cell_type),
          "' and condition '", condition, "'.", call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

.rc_complete_pando_condition_fits_in_object <- function(grn_object) {
  grn <- methods::slot(grn_object, "grn")
  params <- methods::slot(grn, "params")
  fits <- params$condition_grn_fits
  if (!is.list(fits) || !length(fits)) return(grn_object)

  params$condition_grn_fits <-
    .rc_complete_pando_condition_fit_contracts(fits)
  methods::slot(grn, "params") <- params
  methods::slot(grn_object, "grn") <- grn
  grn_object
}

.rc_extract_condition_grn_contract_impl <-
  .rc_extract_condition_grn_contract

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col) {
  grn_object <- .rc_complete_pando_condition_fits_in_object(grn_object)
  grn <- methods::slot(grn_object, "grn")
  params <- methods::slot(grn, "params")
  fits <- params$condition_grn_fits
  fit_list <- if (inherits(fits, "ConditionGRNFit")) list(fits) else fits
  if (!is.list(fit_list) || !length(fit_list)) {
    stop("Pando did not store condition fits for RegCompass extraction.",
         call. = FALSE)
  }
  invisible(lapply(fit_list, .rc_require_pando_condition_grn_fit))

  data_object <- methods::slot(grn_object, "data")
  metadata <- methods::slot(data_object, "meta.data")
  .rc_validate_pando_fit_metadata_frame(
    metadata = metadata,
    fits = fit_list,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  .rc_extract_condition_grn_contract_impl(
    grn_object = grn_object,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
}
