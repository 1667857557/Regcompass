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
  if (!is.list(cells_by_condition) || is.null(names(cells_by_condition)) ||
      !all(levels %in% names(cells_by_condition))) {
    stop("Pando condition cell IDs are not named for every condition.",
         call. = FALSE)
  }
  cells_by_condition <- cells_by_condition[levels]
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
  .rc_extract_condition_grn_contract_impl(
    grn_object = .rc_complete_pando_condition_fits_in_object(grn_object),
    condition_col = condition_col,
    celltype_col = celltype_col
  )
}
