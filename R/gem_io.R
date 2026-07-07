#' Construct a minimal GEM object from stoichiometry and bounds
#' @export
rc_make_gem <- function(S, lb = NULL, ub = NULL, reaction_meta = NULL, metabolite_meta = NULL, medium_policy = "base_bounds", model_info = NULL) {
  S <- as.matrix(S)
  if (is.null(colnames(S))) stop("`S` must have reaction IDs in colnames.", call. = FALSE)
  if (is.null(rownames(S))) rownames(S) <- paste0("met_", seq_len(nrow(S)))
  gem <- list(S = S, lb = lb, ub = ub, reaction_meta = reaction_meta, metabolite_meta = metabolite_meta, medium_policy = medium_policy, model_info = model_info)
  gv <- rc_validate_gem(gem)
  gem$lb <- gv$lb; gem$ub <- gv$ub; gem
}

#' Read a GEM object stored as an RDS file
#' @export
rc_read_gem <- function(file, require_model_info = TRUE) {
  gem <- readRDS(file)
  rc_validate_gem(gem)
  if (isTRUE(require_model_info)) rc_validate_model_info(gem$model_info)
  gem
}

rc_validate_model_info <- function(model_info) {
  if (is.null(model_info) || !is.list(model_info)) {
    stop("GEM is missing `model_info`; set `require_model_info = FALSE` only for legacy exploratory inputs.", call. = FALSE)
  }
  required_legacy <- c("source_model", "model_name", "model_version", "source_release",
                       "source_commit", "file_sha256", "importer", "converted_by",
                       "converted_date")
  required_human2 <- c("source", "version", "commit", "checksum", "conversion_date")
  complete <- function(req) all(vapply(req, function(x) !is.null(model_info[[x]]) && length(model_info[[x]]) > 0L && !all(is.na(model_info[[x]])), logical(1)))
  if (!complete(required_legacy) && !complete(required_human2)) {
    miss <- required_legacy[!vapply(required_legacy, function(x) !is.null(model_info[[x]]) && length(model_info[[x]]) > 0L && !all(is.na(model_info[[x]])), logical(1))]
    stop("`model_info` is incomplete: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}
