# Row-ID parser for labeled microCOMPASS IDs.

rc_parse_microcompass_row_id <- function(x) {
  x <- as.character(x)
  required <- c("reaction", "direction", "medium")

  parse_one <- function(id) {
    fields <- strsplit(id, "::", fixed = TRUE)[[1L]]
    equals <- regexpr("=", fields, fixed = TRUE)
    if (any(equals < 2L)) {
      stop(
        "microCOMPASS row IDs must use `reaction=...::direction=...::medium=...`.",
        call. = FALSE
      )
    }
    keys <- substring(fields, 1L, equals - 1L)
    values <- utils::URLdecode(substring(fields, equals + 1L))
    required_counts <- table(factor(keys, levels = required))
    required_values <- values[match(required, keys)]
    invalid <- any(required_counts != 1L) ||
      anyNA(required_values) ||
      any(!nzchar(trimws(required_values))) ||
      !required_values[[2L]] %in% c("forward", "reverse")
    if (invalid) {
      stop(
        paste(
          "microCOMPASS row IDs must contain exactly one non-empty",
          "`reaction`, `direction`, and `medium` field; `direction` must be",
          "`forward` or `reverse`."
        ),
        call. = FALSE
      )
    }
    named <- stats::setNames(values, keys)
    value <- function(name) {
      hit <- named[name]
      if (!length(hit) || is.na(hit[[1L]]) || !nzchar(trimws(hit[[1L]]))) {
        return(NA_character_)
      }
      unname(hit[[1L]])
    }
    data.frame(
      sample_id = value("sample"),
      module_id = value("module"),
      reaction_id = value("reaction"),
      target_direction = value("direction"),
      medium_scenario = value("medium"),
      condition = value("condition"),
      stringsAsFactors = FALSE
    )
  }

  rows <- lapply(x, parse_one)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
