# Backward compatibility for Pando fits created before the top-level
# dictionary provenance flag was copied from the frozen edge dictionary.

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
