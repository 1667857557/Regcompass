# Exact exchange-to-metabolite correspondence ----------------------------------

.rc_medium_resolve_exchange_metabolites <- function(
    gem, exchange_meta, validated) {
  exchange_ids <- as.character(exchange_meta$reaction_id)
  output <- data.frame(
    exchange_reaction_id = exchange_ids,
    metabolite_id = NA_character_,
    gem_metabolite_name = NA_character_,
    mapping_source = NA_character_,
    stringsAsFactors = FALSE
  )
  explicit_id_col <- intersect(
    c("metabolite_id", "exchange_metabolite_id"),
    colnames(exchange_meta)
  )
  explicit_name_col <- intersect(
    c("metabolite_name", "exchange_metabolite_name"),
    colnames(exchange_meta)
  )
  if (length(explicit_id_col)) {
    output$metabolite_id <- as.character(
      exchange_meta[[explicit_id_col[[1L]]]]
    )
  }
  if (length(explicit_name_col)) {
    output$gem_metabolite_name <- as.character(
      exchange_meta[[explicit_name_col[[1L]]]]
    )
  }
  explicit_id <- !is.na(output$metabolite_id) & nzchar(output$metabolite_id)
  explicit_name <- !is.na(output$gem_metabolite_name) &
    nzchar(output$gem_metabolite_name)
  output$mapping_source[explicit_id | explicit_name] <- "reaction_metadata"

  metabolite_meta <- gem$metabolite_meta
  if (is.null(metabolite_meta) || !is.data.frame(metabolite_meta)) {
    metabolite_meta <- data.frame(
      metabolite_id = validated$metabolites,
      stringsAsFactors = FALSE
    )
  }
  if (!"metabolite_id" %in% colnames(metabolite_meta)) {
    metabolite_meta$metabolite_id <- validated$metabolites
  }
  metabolite_meta <- metabolite_meta[
    match(
      validated$metabolites,
      as.character(metabolite_meta$metabolite_id)
    ),
    ,
    drop = FALSE
  ]
  metabolite_meta$metabolite_id <- validated$metabolites
  name_col <- intersect(
    c("name", "metabolite_name"),
    colnames(metabolite_meta)
  )
  compartment <- if ("compartment" %in% colnames(metabolite_meta)) {
    tolower(as.character(metabolite_meta$compartment))
  } else {
    rep(NA_character_, nrow(metabolite_meta))
  }
  inferred_compartment <- ifelse(
    grepl("\\[e\\]$", validated$metabolites, ignore.case = TRUE),
    "e",
    ifelse(
      grepl("e$", validated$metabolites, ignore.case = TRUE),
      "e",
      NA_character_
    )
  )
  missing_compartment <- is.na(compartment) | !nzchar(compartment)
  compartment[missing_compartment] <-
    inferred_compartment[missing_compartment]

  incomplete <- which(!explicit_id | !explicit_name)
  for (row in incomplete) {
    reaction_index <- match(exchange_ids[[row]], validated$reactions)
    if (is.na(reaction_index)) next
    column <- validated$S[, reaction_index, drop = FALSE]
    nonzero <- which(as.numeric(column) != 0)
    if (!length(nonzero)) next
    extracellular <- nonzero[
      !is.na(compartment[nonzero]) & compartment[nonzero] == "e"
    ]
    candidate <- if (length(extracellular) == 1L) {
      extracellular
    } else if (!length(extracellular) && length(nonzero) == 1L) {
      nonzero
    } else {
      integer()
    }
    if (length(candidate) != 1L) next
    if (!explicit_id[[row]]) {
      output$metabolite_id[[row]] <- validated$metabolites[[candidate]]
    }
    if (!explicit_name[[row]] && length(name_col)) {
      output$gem_metabolite_name[[row]] <- as.character(
        metabolite_meta[[name_col[[1L]]]][[candidate]]
      )
    }
    stoichiometric_source <- if (length(extracellular) == 1L) {
      "stoichiometric_extracellular_metabolite"
    } else {
      "single_stoichiometric_metabolite"
    }
    output$mapping_source[[row]] <- if (
      identical(output$mapping_source[[row]], "reaction_metadata")
    ) {
      paste("reaction_metadata", stoichiometric_source, sep = "+")
    } else {
      stoichiometric_source
    }
  }
  output$normalized_metabolite_id <-
    .rc_medium_normalize_name(output$metabolite_id)
  output$normalized_metabolite_name <-
    .rc_medium_normalize_name(output$gem_metabolite_name)
  output
}

.rc_medium_exchange_metabolites <- .rc_medium_resolve_exchange_metabolites
