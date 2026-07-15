#' Construct a minimal GEM object from stoichiometry and bounds
#' @export
rc_make_gem <- function(S, lb = NULL, ub = NULL,
                        reaction_meta = NULL,
                        metabolite_meta = NULL,
                        medium_policy = "base_bounds",
                        model_info = NULL) {
  S <- .rc_as_dgCMatrix(S)
  if (is.null(colnames(S))) {
    stop("`S` must have reaction IDs in colnames.", call. = FALSE)
  }
  if (is.null(rownames(S))) {
    rownames(S) <- paste0("met_", seq_len(nrow(S)))
  }
  gem <- list(
    S = S,
    lb = lb,
    ub = ub,
    reaction_meta = reaction_meta,
    metabolite_meta = metabolite_meta,
    medium_policy = medium_policy,
    model_info = model_info
  )
  validated <- rc_validate_gem(gem)
  gem$S <- validated$S
  gem$lb <- validated$lb
  gem$ub <- validated$ub
  gem
}

#' Read a GEM object stored as an RDS file
#' @export
rc_read_gem <- function(file) {
  gem <- readRDS(file)
  rc_validate_gem(gem)
  rc_validate_model_info(gem$model_info)
  gem
}

rc_validate_model_info <- function(model_info) {
  if (is.null(model_info) || !is.list(model_info)) {
    stop(
      paste(
        "GEM is missing `model_info`."
      ),
      call. = FALSE
    )
  }
  required_human2 <- c(
    "source", "version", "commit", "checksum",
    "conversion_date"
  )
  complete <- function(required) {
    all(vapply(required, function(name) {
      !is.null(model_info[[name]]) &&
        length(model_info[[name]]) > 0L &&
        !all(is.na(model_info[[name]]))
    }, logical(1)))
  }
  if (!complete(required_human2)) {
    missing <- required_human2[!vapply(required_human2, function(name) {
      !is.null(model_info[[name]]) &&
        length(model_info[[name]]) > 0L &&
        !all(is.na(model_info[[name]]))
    }, logical(1))]
    stop(
      "`model_info` is incomplete: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
