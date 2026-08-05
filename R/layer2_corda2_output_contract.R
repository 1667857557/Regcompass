# Preserve both directional-variable and reaction-level CORDA2 associations.

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
