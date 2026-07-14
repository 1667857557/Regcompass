# Row-ID parser supporting both v1.3 labeled IDs and legacy IDs.

rc_parse_microcompass_row_id <- function(x) {
  x <- as.character(x)
  labeled <- grepl("^sample=", x)
  output <- data.frame(
    sample_id = rep(NA_character_, length(x)),
    module_id = rep(NA_character_, length(x)),
    reaction_id = rep(NA_character_, length(x)),
    target_direction = rep(NA_character_, length(x)),
    medium_scenario = rep(NA_character_, length(x)),
    stringsAsFactors = FALSE
  )

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

  if (any(labeled)) {
    parsed <- lapply(x[labeled], function(value) {
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
    output$sample_id[labeled] <- vapply(
      parsed, get_named_value, character(1), "sample"
    )
    output$module_id[labeled] <- vapply(
      parsed, get_named_value, character(1), "module"
    )
    output$reaction_id[labeled] <- vapply(
      parsed, get_named_value, character(1), "reaction"
    )
    output$target_direction[labeled] <- vapply(
      parsed, get_named_value, character(1), "direction"
    )
    output$medium_scenario[labeled] <- vapply(
      parsed, get_named_value, character(1), "medium"
    )
  }

  if (any(!labeled)) {
    legacy <- x[!labeled]
    core <- sub("::medium=.*$", "", legacy)
    parts <- strsplit(core, "::", fixed = TRUE)
    output$reaction_id[!labeled] <- vapply(
      parts, get_position, character(1), 1L
    )
    output$target_direction[!labeled] <- vapply(
      parts, get_position, character(1), 2L
    )
    output$medium_scenario[!labeled] <- ifelse(
      grepl("::medium=", legacy, fixed = TRUE),
      sub("^.*::medium=", "", legacy),
      vapply(parts, get_position, character(1), 3L)
    )
  }
  output
}
