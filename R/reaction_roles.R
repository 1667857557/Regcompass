.rc_meta_col <- function(meta, names, n) {
  column <- intersect(names, colnames(meta))[1L]
  if (is.na(column)) return(rep("", n))
  value <- as.character(meta[[column]])
  value[is.na(value)] <- ""
  tolower(value)
}

.rc_infer_exchange_like <- function(gem, validated, meta, role,
                                    medium_table = NULL,
                                    infer_from_id = TRUE,
                                    infer_from_stoichiometry = TRUE,
                                    infer_from_compartment = TRUE) {
  n <- length(validated$reactions)
  evidence <- rep(FALSE, n)
  reaction_id <- tolower(as.character(validated$reactions))

  if (isTRUE(infer_from_id)) {
    reaction_name <- .rc_meta_col(
      meta, c("reaction_name", "name", "description"), n
    )
    subsystem <- .rc_meta_col(
      meta,
      c("subsystem", "subSystems", "subSystem", "sub_system", "metabolic_module"),
      n
    )
    evidence <- evidence |
      grepl("^(ex_|ex-|exchange)", reaction_id) |
      grepl("exchange|boundary exchange|extracellular exchange|uptake|secretion", reaction_name) |
      grepl("exchange|uptake|secretion", subsystem)
  }

  if (isTRUE(infer_from_stoichiometry)) {
    equation <- .rc_meta_col(
      meta, c("equation", "reaction_formula", "formula"), n
    )
    boundary_like <- Matrix::colSums(abs(validated$S) > 0) == 1L
    extracellular_equation <- grepl(
      "\\[e\\]|\\[extracellular\\]|extracellular", equation
    )
    evidence <- evidence | (boundary_like & extracellular_equation)
  }

  if (isTRUE(infer_from_compartment) &&
      !is.null(gem$metabolite_meta) &&
      all(c("metabolite_id", "compartment") %in% colnames(gem$metabolite_meta))) {
    compartment <- as.character(gem$metabolite_meta$compartment[
      match(
        rownames(validated$S),
        as.character(gem$metabolite_meta$metabolite_id)
      )
    ])
    single_extracellular <- vapply(seq_len(ncol(validated$S)), function(j) {
      index <- which(validated$S[, j] != 0)
      reaction_compartment <- unique(stats::na.omit(compartment[index]))
      length(index) == 1L && length(reaction_compartment) == 1L &&
        identical(tolower(reaction_compartment[[1L]]), "e")
    }, logical(1))
    evidence <- evidence | single_extracellular
  }

  if (!is.null(medium_table) &&
      is.data.frame(medium_table) &&
      "exchange_reaction_id" %in% colnames(medium_table)) {
    evidence <- evidence |
      validated$reactions %in% as.character(medium_table$exchange_reaction_id)
  }
  role == "unknown" & evidence
}

#' Annotate GEM reactions with curated or inferred RegCompassR roles
#' @export
rc_annotate_reaction_roles <- function(
    gem, reaction_role_table = NULL, medium_table = NULL,
    infer_from_id = TRUE, infer_from_stoichiometry = TRUE,
    infer_from_compartment = TRUE, overwrite_existing = FALSE) {
  flags <- c(
    infer_from_id = infer_from_id,
    infer_from_stoichiometry = infer_from_stoichiometry,
    infer_from_compartment = infer_from_compartment,
    overwrite_existing = overwrite_existing
  )
  if (any(vapply(flags, function(x) {
    !is.logical(x) || length(x) != 1L || is.na(x)
  }, logical(1)))) {
    stop("Reaction-role switches must be TRUE or FALSE.", call. = FALSE)
  }

  validated <- rc_validate_gem(gem)
  meta <- gem$reaction_meta
  if (is.null(meta) || !is.data.frame(meta)) {
    meta <- data.frame(
      reaction_id = validated$reactions,
      stringsAsFactors = FALSE
    )
  }
  if (!"reaction_id" %in% colnames(meta)) {
    if (nrow(meta) != length(validated$reactions)) {
      stop(
        "Reaction metadata without `reaction_id` must have one row per GEM reaction.",
        call. = FALSE
      )
    }
    meta$reaction_id <- validated$reactions
  }
  metadata_ids <- trimws(as.character(meta$reaction_id))
  if (anyNA(metadata_ids) || any(!nzchar(metadata_ids)) ||
      anyDuplicated(metadata_ids)) {
    stop("Reaction metadata IDs must be unique and non-empty.", call. = FALSE)
  }
  missing_meta <- setdiff(validated$reactions, metadata_ids)
  if (length(missing_meta)) {
    filler <- as.data.frame(
      matrix(NA, nrow = length(missing_meta), ncol = ncol(meta)),
      stringsAsFactors = FALSE
    )
    names(filler) <- names(meta)
    filler$reaction_id <- missing_meta
    meta <- rbind(meta, filler)
    metadata_ids <- as.character(meta$reaction_id)
  }
  meta <- meta[
    match(validated$reactions, metadata_ids),
    , drop = FALSE
  ]
  meta$reaction_id <- validated$reactions

  old_role <- if ("role" %in% colnames(meta)) {
    trimws(as.character(meta$role))
  } else {
    rep(NA_character_, nrow(meta))
  }
  role <- ifelse(
    !is.na(old_role) & nzchar(old_role) & !overwrite_existing,
    old_role,
    "unknown"
  )
  source <- ifelse(role != "unknown", "metadata", "unknown")
  confidence <- ifelse(role != "unknown", "medium", "low")

  exchange_index <- .rc_infer_exchange_like(
    gem = gem,
    validated = validated,
    meta = meta,
    role = role,
    medium_table = medium_table,
    infer_from_id = infer_from_id,
    infer_from_stoichiometry = infer_from_stoichiometry,
    infer_from_compartment = infer_from_compartment
  )
  role[exchange_index] <- "exchange"
  source[exchange_index] <- "inferred_exchange_evidence"
  confidence[exchange_index] <- "medium"

  if (infer_from_stoichiometry) {
    nonzero <- Matrix::colSums(abs(validated$S) > 0)
    index <- role == "unknown" & nonzero == 1L
    role[index] <- "boundary_like"
    source[index] <- "stoichiometry"
    confidence[index] <- "low"
  }

  if (infer_from_compartment &&
      !is.null(gem$metabolite_meta) &&
      all(c("metabolite_id", "compartment") %in% colnames(gem$metabolite_meta))) {
    compartment <- as.character(gem$metabolite_meta$compartment[
      match(
        rownames(validated$S),
        as.character(gem$metabolite_meta$metabolite_id)
      )
    ])
    compartments_by_reaction <- lapply(seq_len(ncol(validated$S)), function(j) {
      unique(stats::na.omit(compartment[which(validated$S[, j] != 0)]))
    })
    index <- role == "unknown" &
      vapply(compartments_by_reaction, length, integer(1)) > 1L
    role[index] <- "transport"
    source[index] <- "compartment_stoichiometry"
    confidence[index] <- "medium"
  }

  if (infer_from_id) {
    reaction_id <- tolower(validated$reactions)
    patterns <- list(
      exchange = "^(ex_|ex-|exchange)",
      demand = "^(dm_|demand)",
      sink = "^(sink_|sk_)",
      biomass = "biomass",
      maintenance = "maintenance|atpm",
      transport = "(_tx$|_transport$|transport)"
    )
    for (name in names(patterns)) {
      index <- role == "unknown" & grepl(patterns[[name]], reaction_id)
      role[index] <- name
      source[index] <- "id_pattern"
      confidence[index] <- "low"
    }
  }

  if (!is.null(reaction_role_table)) {
    if (!is.data.frame(reaction_role_table) ||
        !all(c("reaction_id", "role") %in% colnames(reaction_role_table))) {
      stop(
        "`reaction_role_table` must contain `reaction_id` and `role`.",
        call. = FALSE
      )
    }
    curated_id <- trimws(as.character(reaction_role_table$reaction_id))
    curated_role <- trimws(as.character(reaction_role_table$role))
    if (anyNA(curated_id) || any(!nzchar(curated_id)) ||
        anyDuplicated(curated_id) || anyNA(curated_role) ||
        any(!nzchar(curated_role))) {
      stop("Curated reaction roles require unique, non-empty IDs and roles.", call. = FALSE)
    }
    unknown <- setdiff(curated_id, validated$reactions)
    if (length(unknown)) {
      stop(
        "Curated reaction roles contain unknown reactions: ",
        paste(utils::head(unknown, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    match_index <- match(validated$reactions, curated_id)
    index <- !is.na(match_index)
    role[index] <- curated_role[match_index[index]]
    source[index] <- "curated"
    confidence[index] <- "high"
  }

  meta$role <- role
  meta$role_source <- source
  meta$role_confidence <- confidence
  gem$reaction_meta <- meta
  gem$reaction_roles <- meta[
    , c("reaction_id", "role", "role_source", "role_confidence"),
    drop = FALSE
  ]
  gem
}
