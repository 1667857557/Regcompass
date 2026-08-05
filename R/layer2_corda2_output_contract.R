# Preserve directional CORDA2 associations and validate early metadata.

.rc_validate_corda_union_model_before_corda2_output <-
  .rc_validate_corda_union_model

.rc_validate_corda_union_model <- function(model, cell_type) {
  algorithm <- as.character(model$build_params$algorithm %||% "")
  if (identical(
    algorithm,
    "resendislab_python_CORDA2_corrected_redundant_path_assessment"
  )) {
    reference <-
      model$build_params$python_reference_commit %||%
      model$corda_reconstruction$python_reference_commit %||%
      "c02e06d50606bf93f23d8f2e6d6ade0e996ca70e"
    model$build_params$python_reference_commit <- reference
  }
  .rc_validate_corda_union_model_before_corda2_output(
    model, cell_type
  )
}

.rc_corda_normalize_associations_before_corda2 <-
  .rc_corda_normalize_associations

.rc_corda2_reaction_from_variable <- function(variable_id) {
  sub("::(forward|reverse)$", "", as.character(variable_id))
}

.rc_corda_normalize_associations <- function(tab) {
  answer <- .rc_corda_normalize_associations_before_corda2(tab)
  if (!is.data.frame(answer)) return(answer)
  if (!nrow(answer)) {
    answer$associated_variable_id <- character()
    return(answer[, c(
      "task_key", "associated_variable_id", "associated_reaction_id"
    ), drop = FALSE])
  }
  variable <- as.character(answer$associated_reaction_id)
  answer$associated_variable_id <- variable
  answer$associated_reaction_id <-
    .rc_corda2_reaction_from_variable(variable)
  answer <- answer[, c(
    "task_key", "associated_variable_id", "associated_reaction_id"
  ), drop = FALSE]
  answer[order(
    answer$task_key,
    answer$associated_reaction_id,
    answer$associated_variable_id
  ), , drop = FALSE]
}

.rc_complete_celltype_medium_corda_gem_before_corda2_output <-
  .rc_complete_celltype_medium_corda_gem

.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  model <- do.call(
    .rc_complete_celltype_medium_corda_gem_before_corda2_output,
    args
  )
  algorithm <- as.character(model$build_params$algorithm %||% "")
  if (!identical(
    algorithm,
    "resendislab_python_CORDA2_corrected_redundant_path_assessment"
  )) {
    return(model)
  }
  edges <- model$corda_association_edges
  if (!is.data.frame(edges)) {
    edges <- data.frame(
      task_key = character(),
      associated_variable_id = character(),
      associated_reaction_id = character(),
      stringsAsFactors = FALSE
    )
  } else {
    if (!"associated_variable_id" %in% colnames(edges)) {
      edges$associated_variable_id <- if (nrow(edges)) {
        as.character(edges$associated_reaction_id)
      } else {
        character()
      }
    }
    edges <- edges[, c(
      "task_key", "associated_variable_id", "associated_reaction_id"
    ), drop = FALSE]
  }
  model$corda_association_edges <- edges
  model$build_params$association_edge_schema <- c(
    "task_key", "associated_variable_id", "associated_reaction_id"
  )
  .rc_validate_corda_union_model(
    model,
    cell_type = as.character(args$cell_type)
  )
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
