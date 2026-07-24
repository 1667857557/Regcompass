.rc_stage_gem_fingerprint <- function(gem) {
  .rc_full_gem_cache_fingerprint(gem)
}

.rc_require_stage_class <- function(x, class_name, argument, producer) {
  if (!inherits(x, class_name)) {
    stop(
      "`", argument, "` must be the output of `", producer, "()`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_require_stage_gem <- function(x, gem, argument) {
  expected <- as.character(x$gem_fingerprint %||% "")
  observed <- .rc_stage_gem_fingerprint(gem)
  if (length(expected) != 1L || !nzchar(expected)) {
    stop("`", argument, "` lacks GEM provenance.", call. = FALSE)
  }
  if (!identical(expected, observed)) {
    stop("`", argument, "` was generated from a different GEM.", call. = FALSE)
  }
  invisible(observed)
}

.rc_require_workflow_params <- function(x, expected, argument) {
  observed <- x$workflow_params %||% x$params
  if (!is.list(observed) || !identical(observed, expected)) {
    stop("`", argument, "` uses different workflow parameters.", call. = FALSE)
  }
  invisible(TRUE)
}

.rc_layer1_unit_ids <- function(layer1) {
  expression <- layer1$reaction_expression
  meta <- layer1$unit_meta
  if (!is.numeric(expression) || is.null(dim(expression)) ||
      is.null(rownames(expression)) || is.null(colnames(expression)) ||
      anyNA(rownames(expression)) || anyNA(colnames(expression)) ||
      any(!nzchar(rownames(expression))) || any(!nzchar(colnames(expression))) ||
      anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression))) {
    stop(
      "Layer 1 reaction expression must be a numeric matrix with unique reaction and unit IDs.",
      call. = FALSE
    )
  }
  if (!is.data.frame(meta)) {
    stop("Layer 1 `unit_meta` must be a data frame.", call. = FALSE)
  }
  id_col <- if ("pool_id" %in% colnames(meta)) {
    "pool_id"
  } else if ("unit_id" %in% colnames(meta)) {
    "unit_id"
  } else {
    stop("Layer 1 `unit_meta` lacks `pool_id`/`unit_id`.", call. = FALSE)
  }
  ids <- trimws(as.character(meta[[id_col]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Layer 1 unit IDs must be unique and non-empty.", call. = FALSE)
  }
  if (!identical(colnames(expression), ids)) {
    stop(
      "Layer 1 reaction-expression columns and `unit_meta` are not identically ordered.",
      call. = FALSE
    )
  }
  ids
}

.rc_validate_layer1_stage <- function(
    layer1, workflow_params = NULL, gem = NULL,
    argument = "layer1") {
  .rc_require_stage_class(
    layer1, "regcompass_layer1_step", argument,
    "rc_regcompass_step_layer1"
  )
  ids <- .rc_layer1_unit_ids(layer1)
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer1, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer1, gem, argument)
  invisible(ids)
}

.rc_layer2_unit_ids <- function(layer2) {
  required <- c("penalty", "vmax", "feasible", "evaluated", "score")
  missing <- required[!vapply(required, function(name) {
    !is.null(layer2[[name]]) && !is.null(dim(layer2[[name]]))
  }, logical(1))]
  if (length(missing)) {
    stop("Layer 2 is missing matrices: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  reference <- layer2$penalty
  if (!is.numeric(reference) || is.null(rownames(reference)) ||
      is.null(colnames(reference)) || anyDuplicated(rownames(reference)) ||
      anyDuplicated(colnames(reference))) {
    stop("Layer 2 penalty requires unique target and unit IDs.", call. = FALSE)
  }
  for (name in setdiff(required, "penalty")) {
    value <- layer2[[name]]
    if (!identical(dimnames(value), dimnames(reference))) {
      stop("Layer 2 `", name, "` is not aligned with `penalty`.", call. = FALSE)
    }
  }
  meta <- layer2$unit_meta
  if (!is.data.frame(meta)) {
    stop("Layer 2 `unit_meta` must be a data frame.", call. = FALSE)
  }
  id_col <- if ("pool_id" %in% colnames(meta)) {
    "pool_id"
  } else if ("unit_id" %in% colnames(meta)) {
    "unit_id"
  } else {
    stop("Layer 2 `unit_meta` lacks `pool_id`/`unit_id`.", call. = FALSE)
  }
  ids <- trimws(as.character(meta[[id_col]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      !identical(colnames(reference), ids)) {
    stop("Layer 2 matrices and `unit_meta` are not identically aligned.",
         call. = FALSE)
  }
  ids
}

.rc_validate_layer2_stage <- function(
    layer2, layer1 = NULL, workflow_params = NULL, gem = NULL,
    required_mode = NULL, argument = "layer2") {
  .rc_require_stage_class(
    layer2, "regcompass_layer2_step", argument,
    "rc_regcompass_step_layer2"
  )
  ids <- .rc_layer2_unit_ids(layer2)
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer2, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer2, gem, argument)
  if (!is.null(required_mode) &&
      !identical(as.character(layer2$model_mode), required_mode)) {
    stop("`", argument, "` must use `model_mode = \"", required_mode, "\"`.",
         call. = FALSE)
  }
  if (!is.null(layer1)) {
    layer1_ids <- .rc_validate_layer1_stage(
      layer1,
      workflow_params = workflow_params,
      gem = gem,
      argument = "layer1"
    )
    if (!identical(ids, layer1_ids)) {
      stop("Layer 1 and Layer 2 contain different scoring units.", call. = FALSE)
    }
  }
  invisible(ids)
}
