# Row-ID parser for v1.3 labeled microCOMPASS IDs.

rc_parse_microcompass_row_id <- function(x) {
  x <- as.character(x)
  required_fields <- c("reaction", "direction", "medium")

  get_named_value <- function(value, name) {
    selected <- value[name]
    if (!length(selected) || is.na(selected[[1L]]) ||
        !nzchar(selected[[1L]])) {
      return(NA_character_)
    }
    unname(selected[[1L]])
  }
  get_position <- function(value, position) {
    if (length(value) < position) return(NA_character_)
    value[[position]]
  }

  has_required_labels <- vapply(x, function(value) {
    fields <- strsplit(value, "::", fixed = TRUE)[[1L]]
    keys <- vapply(
      strsplit(fields, "=", fixed = TRUE),
      get_position,
      character(1),
      1L
    )
    all(required_fields %in% keys)
  }, logical(1))
  if (any(!has_required_labels)) {
    stop(
      paste(
        "microCOMPASS row IDs must use the v1.3 labeled format",
        "`reaction=...::direction=...::medium=...`."
      ),
      call. = FALSE
    )
  }
  output <- data.frame(
    sample_id = rep(NA_character_, length(x)),
    module_id = rep(NA_character_, length(x)),
    reaction_id = rep(NA_character_, length(x)),
    target_direction = rep(NA_character_, length(x)),
    medium_scenario = rep(NA_character_, length(x)),
    condition = rep(NA_character_, length(x)),
    stringsAsFactors = FALSE
  )

  parsed <- lapply(x, function(value) {
    fields <- strsplit(value, "::", fixed = TRUE)[[1L]]
    key_value <- strsplit(fields, "=", fixed = TRUE)
    keys <- vapply(key_value, get_position, character(1), 1L)
    values <- vapply(key_value, function(field) {
      if (length(field) < 2L) return(NA_character_)
      paste(field[-1L], collapse = "=")
    }, character(1))
    values <- utils::URLdecode(values)
    names(values) <- keys
    values
  })
  invalid <- vapply(parsed, function(values) {
    duplicated_required <- any(duplicated(names(values)[names(values) %in% required_fields]))
    missing_required <- any(vapply(required_fields, function(name) {
      value <- get_named_value(values, name)
      is.na(value) || !nzchar(value)
    }, logical(1)))
    direction <- get_named_value(values, "direction")
    bad_direction <- is.na(direction) || !direction %in% c("forward", "reverse")
    duplicated_required || missing_required || bad_direction
  }, logical(1))
  if (any(invalid)) {
    stop(
      paste(
        "microCOMPASS row IDs must contain exactly one non-empty",
        "`reaction`, `direction`, and `medium` field; `direction` must be",
        "`forward` or `reverse`."
      ),
      call. = FALSE
    )
  }
  output$sample_id <- vapply(parsed, get_named_value, character(1), "sample")
  output$module_id <- vapply(parsed, get_named_value, character(1), "module")
  output$reaction_id <- vapply(parsed, get_named_value, character(1), "reaction")
  output$target_direction <- vapply(parsed, get_named_value, character(1), "direction")
  output$medium_scenario <- vapply(parsed, get_named_value, character(1), "medium")
  output$condition <- vapply(parsed, get_named_value, character(1), "condition")
  output
}
