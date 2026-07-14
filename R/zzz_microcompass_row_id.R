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

  if (any(labeled)) {
    parsed <- lapply(x[labeled], function(value) {
      fields <- strsplit(value, "::", fixed = TRUE)[[1L]]
      key_value <- strsplit(fields, "=", fixed = TRUE)
      keys <- vapply(key_value, `[[`, character(1), 1L)
      values <- vapply(key_value, function(field) {
        if (length(field) < 2L) return(NA_character_)
        paste(field[-1L], collapse = "=")
      }, character(1))
      values <- utils::URLdecode(values)
      names(values) <- keys
      values
    })
    output$sample_id[labeled] <- vapply(
      parsed,
      function(value) value[["sample"]] %||% NA_character_,
      character(1)
    )
    output$module_id[labeled] <- vapply(
      parsed,
      function(value) value[["module"]] %||% NA_character_,
      character(1)
    )
    output$reaction_id[labeled] <- vapply(
      parsed,
      function(value) value[["reaction"]] %||% NA_character_,
      character(1)
    )
    output$target_direction[labeled] <- vapply(
      parsed,
      function(value) value[["direction"]] %||% NA_character_,
      character(1)
    )
    output$medium_scenario[labeled] <- vapply(
      parsed,
      function(value) value[["medium"]] %||% NA_character_,
      character(1)
    )
  }

  if (any(!labeled)) {
    legacy <- x[!labeled]
    core <- sub("::medium=.*$", "", legacy)
    parts <- strsplit(core, "::", fixed = TRUE)
    output$reaction_id[!labeled] <- vapply(
      parts,
      function(value) value[[1L]] %||% NA_character_,
      character(1)
    )
    output$target_direction[!labeled] <- vapply(
      parts,
      function(value) value[[2L]] %||% NA_character_,
      character(1)
    )
    output$medium_scenario[!labeled] <- ifelse(
      grepl("::medium=", legacy, fixed = TRUE),
      sub("^.*::medium=", "", legacy),
      vapply(
        parts,
        function(value) value[[3L]] %||% NA_character_,
        character(1)
      )
    )
  }
  output
}
