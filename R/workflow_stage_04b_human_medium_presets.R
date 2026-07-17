# Workflow stage 4b: published human medium presets.
#
# This file extends the shared model-bound medium contract with explicit human
# plasma, RPMI-1640, glucose/lactate sensitivity, and custom-metabolite presets.
# Unlisted exchanges remain closed for uptake when the medium is applied.

.rc_medium_presets_previous_make_medium_scenarios <- rc_make_medium_scenarios

.rc_medium_reference_catalog <- function() {
  data.frame(
    preset_id = c(
      "normal_human_plasma", "high_glucose", "low_glucose",
      "high_lactate", "low_lactate", "rpmi1640", "custom"
    ),
    reference_label = c(
      "Cantor et al., Cell 2017 (HPLM); Psychogios et al., PLoS One 2011",
      "Han et al., Gynecologic Oncology 2015",
      "Han et al., Gynecologic Oncology 2015",
      "Schwickert et al., Experientia 1996; Kennedy et al., PLoS One 2013",
      "Kennedy et al., PLoS One 2013",
      "Moore et al., JAMA 1967 (RPMI-1640); Cantor et al., Cell 2017",
      "user supplied"
    ),
    reference_doi = c(
      "10.1016/j.cell.2017.03.023;10.1371/journal.pone.0016957",
      "10.1016/j.ygyno.2015.06.036",
      "10.1016/j.ygyno.2015.06.036",
      "10.1007/BF01919316;10.1371/journal.pone.0075154",
      "10.1371/journal.pone.0075154",
      "10.1001/jama.1967.03120080053007;10.1016/j.cell.2017.03.023",
      NA_character_
    ),
    reference_pmid = c(
      "28388410;21359215", "26135947", "26135947",
      "8641383;24069390", "24069390", "28388410", NA_character_
    ),
    species = "Homo sapiens",
    evidence_scope = c(
      "adult human plasma/serum nutrient availability",
      "25 mM glucose sensitivity condition",
      "1 mM glucose sensitivity condition",
      "20 mM lactate sensitivity condition",
      "0.5 mM lactate sensitivity condition",
      "serum-free RPMI-1640 basal formulation",
      "user-supplied metabolite availability"
    ),
    stringsAsFactors = FALSE
  )
}

.rc_medium_pattern <- function(name) {
  patterns <- c(
    glucose = "glucose|(^|[^a-z0-9])glc([^a-z0-9]|$)",
    lactate = "lactate|lactic[ _-]?acid|(^|[^a-z0-9])lac([^a-z0-9]|$)",
    glutamine = "glutamine|(^|[^a-z0-9])gln([^a-z0-9]|$)",
    arginine = "arginine|(^|[^a-z0-9])arg([^a-z0-9]|$)",
    oxygen = "oxygen|(^|[^a-z0-9])o2([^a-z0-9]|$)",
    urate = "urate|uric[ _-]?acid"
  )
  if (name %in% names(patterns)) return(unname(patterns[[name]]))
  gsub("_", "[- _]?", name, fixed = TRUE)
}

.rc_medium_compounds <- function(preset_id) {
  plasma <- data.frame(
    metabolite_name = c(
      "glucose", "lactate", "glutamine", "arginine", "oxygen", "urate"
    ),
    concentration_mM = c(5, 1.5, 0.55, 0.11, NA, 0.35),
    uptake_fraction = c(5 / 25, 1.5 / 20, 1, 1, 1, 1),
    target_exchange_flag = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    required_match = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  plasma$metabolite_pattern <- vapply(
    plasma$metabolite_name, .rc_medium_pattern, character(1)
  )

  if (identical(preset_id, "normal_human_plasma")) return(plasma)
  if (identical(preset_id, "rpmi1640")) {
    rpm <- data.frame(
      metabolite_name = c("glucose", "glutamine", "arginine", "oxygen"),
      concentration_mM = c(11.1, 2.055, 1.149, NA),
      uptake_fraction = c(11.1 / 25, 1, 1, 1),
      target_exchange_flag = c(TRUE, FALSE, FALSE, FALSE),
      required_match = c(TRUE, TRUE, FALSE, FALSE),
      stringsAsFactors = FALSE
    )
    rpm$metabolite_pattern <- vapply(
      rpm$metabolite_name, .rc_medium_pattern, character(1)
    )
    return(rpm)
  }

  target <- if (preset_id %in% c("high_glucose", "low_glucose")) {
    "glucose"
  } else if (preset_id %in% c("high_lactate", "low_lactate")) {
    "lactate"
  } else {
    stop("Unsupported human medium preset: ", preset_id, call. = FALSE)
  }
  concentration <- switch(
    preset_id,
    high_glucose = 25,
    low_glucose = 1,
    high_lactate = 20,
    low_lactate = 0.5
  )
  reference_high <- if (identical(target, "glucose")) 25 else 20
  selected <- plasma$metabolite_name == target
  plasma$concentration_mM[selected] <- concentration
  plasma$uptake_fraction[selected] <- concentration / reference_high
  plasma$target_exchange_flag[selected] <- TRUE
  plasma$required_match[selected] <- TRUE
  plasma
}

.rc_medium_annotation_text <- function(meta) {
  columns <- intersect(
    c(
      "reaction_id", "reaction_name", "name", "description", "equation",
      "metabolite_id", "metabolite_name"
    ),
    colnames(meta)
  )
  if (!length(columns)) return(tolower(as.character(meta$reaction_id)))
  values <- meta[, columns, drop = FALSE]
  values[] <- lapply(values, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  tolower(apply(values, 1L, paste, collapse = " "))
}

.rc_medium_scale <- function(uptake_scale, scenario_id) {
  if (!is.numeric(uptake_scale) || !length(uptake_scale) ||
      any(!is.finite(uptake_scale)) || any(uptake_scale < 0)) {
    stop("`uptake_scale` must contain finite non-negative values.", call. = FALSE)
  }
  if (length(uptake_scale) == 1L &&
      (is.null(names(uptake_scale)) || !nzchar(names(uptake_scale)[[1L]]))) {
    return(as.numeric(uptake_scale[[1L]]))
  }
  if (!is.null(names(uptake_scale)) && scenario_id %in% names(uptake_scale)) {
    return(as.numeric(uptake_scale[[scenario_id]]))
  }
  1
}

.rc_build_human_medium_preset <- function(
    gem, preset_id, exchange_limit, uptake_scale,
    condition = "all", exchange_roles = "exchange",
    strict_preset_matching = TRUE, compounds = NULL,
    custom_reference = NULL) {
  if (!is.numeric(exchange_limit) || length(exchange_limit) != 1L ||
      !is.finite(exchange_limit) || exchange_limit <= 0) {
    stop("`exchange_limit` must be one positive finite number.", call. = FALSE)
  }
  if (!is.logical(strict_preset_matching) ||
      length(strict_preset_matching) != 1L || is.na(strict_preset_matching)) {
    stop("`strict_preset_matching` must be TRUE or FALSE.", call. = FALSE)
  }

  validated <- rc_validate_gem(gem)
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem)
  }
  meta <- gem$reaction_meta[
    match(validated$reactions, as.character(gem$reaction_meta$reaction_id)),
    , drop = FALSE
  ]
  roles <- unique(trimws(as.character(exchange_roles)))
  roles <- roles[!is.na(roles) & nzchar(roles)]
  exchange_meta <- meta[as.character(meta$role) %in% roles, , drop = FALSE]
  if (!nrow(exchange_meta)) {
    stop("No `exchange` reactions found for the medium preset.", call. = FALSE)
  }
  exchange_ids <- as.character(exchange_meta$reaction_id)

  if (is.null(compounds)) compounds <- .rc_medium_compounds(preset_id)
  required <- c(
    "metabolite_name", "metabolite_pattern", "concentration_mM",
    "uptake_fraction", "target_exchange_flag", "required_match"
  )
  missing <- setdiff(required, colnames(compounds))
  if (length(missing)) {
    stop("Preset compound table missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  annotation <- .rc_medium_annotation_text(exchange_meta)
  matched_count <- integer(nrow(compounds))
  rows <- vector("list", nrow(compounds))
  scenario_scale <- .rc_medium_scale(uptake_scale, preset_id)
  for (i in seq_len(nrow(compounds))) {
    hit <- grepl(
      as.character(compounds$metabolite_pattern[[i]]),
      annotation,
      ignore.case = TRUE,
      perl = TRUE
    )
    matched_count[[i]] <- sum(hit)
    if (!any(hit)) next
    index <- match(exchange_ids[hit], validated$reactions)
    fraction <- pmin(
      1,
      as.numeric(compounds$uptake_fraction[[i]]) * scenario_scale
    )
    original_lb <- as.numeric(validated$lb[index])
    original_ub <- as.numeric(validated$ub[index])
    lb <- pmax(original_lb, -exchange_limit * fraction)
    lb[original_lb >= 0] <- 0
    rows[[i]] <- data.frame(
      medium_scenario_id = preset_id,
      exchange_reaction_id = exchange_ids[hit],
      metabolite_id = if ("metabolite_id" %in% colnames(exchange_meta)) {
        as.character(exchange_meta$metabolite_id[hit])
      } else {
        NA_character_
      },
      preset_metabolite = as.character(compounds$metabolite_name[[i]]),
      concentration_mM = as.numeric(compounds$concentration_mM[[i]]),
      condition = as.character(condition %||% "all"),
      lb = lb,
      ub = pmin(original_ub, exchange_limit),
      available = TRUE,
      original_lb = original_lb,
      original_ub = original_ub,
      exchange_limit = exchange_limit,
      uptake_fraction = fraction,
      evidence_source = "published_human_medium_preset",
      assumption_level =
        "literature_defined_availability_with_relative_uptake_cap",
      target_exchange_flag = as.logical(compounds$target_exchange_flag[[i]]),
      concentration_used_for_rate_bound =
        as.logical(compounds$target_exchange_flag[[i]]),
      rate_bound_source = if (isTRUE(compounds$target_exchange_flag[[i]])) {
        "relative_concentration_sensitivity_not_measured_flux"
      } else {
        "binary_preset_availability_not_measured_flux"
      },
      stringsAsFactors = FALSE
    )
  }

  unmatched <- compounds$metabolite_name[
    as.logical(compounds$required_match) & matched_count == 0L
  ]
  if (length(unmatched) && isTRUE(strict_preset_matching)) {
    stop(
      "Required preset metabolites were not matched to GEM exchanges: ",
      paste(unmatched, collapse = ", "),
      ". Inspect `gem$reaction_meta` or use `custom_medium`.",
      call. = FALSE
    )
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    stop("The medium preset did not match any exchange reactions.", call. = FALSE)
  }
  output <- do.call(rbind, rows)
  output <- output[!duplicated(output$exchange_reaction_id), , drop = FALSE]

  reference <- custom_reference %||% .rc_medium_reference_catalog()
  ref <- reference[reference$preset_id == preset_id, , drop = FALSE]
  if (!nrow(ref)) ref <- .rc_medium_reference_catalog()[7L, , drop = FALSE]
  output$reference_label <- ref$reference_label[[1L]]
  output$reference_doi <- ref$reference_doi[[1L]]
  output$reference_pmid <- ref$reference_pmid[[1L]]
  output$species <- ref$species[[1L]]
  output$evidence_scope <- ref$evidence_scope[[1L]]
  rownames(output) <- NULL
  attr(output, "preset_diagnostics") <- data.frame(
    medium_scenario_id = preset_id,
    preset_metabolite = compounds$metabolite_name,
    n_exchange_matches = matched_count,
    required_match = compounds$required_match,
    matched = matched_count > 0L,
    stringsAsFactors = FALSE
  )
  output
}

rc_make_medium_scenarios <- function(
    gem,
    scenario = "compass_model_bounds",
    custom_medium = NULL,
    custom_metabolites = NULL,
    uptake_scale = c(
      permissive_all_exchange = 1,
      normal_human_plasma = 1,
      rpmi1640 = 1,
      minimal = 0.1,
      low_glucose = 0.1,
      low_glutamine = 0.1,
      high_lactate = 1
    ),
    condition_col = NULL,
    exchange_roles = c("exchange"),
    condition = condition_col,
    exchange_limit = 1,
    strict_preset_matching = TRUE) {
  choices <- c(
    "compass_model_bounds", "permissive_all_exchange",
    "normal_human_plasma", "high_glucose", "low_glucose",
    "high_lactate", "low_lactate", "rpmi1640",
    "minimal", "low_glutamine",
    "blood_like", "culture_like", "tumor_low_glucose",
    "lactate_available", "custom"
  )
  scenario <- match.arg(scenario, choices = choices, several.ok = TRUE)

  aliases <- c(
    blood_like = "normal_human_plasma",
    culture_like = "rpmi1640",
    tumor_low_glucose = "low_glucose",
    lactate_available = "high_lactate"
  )
  alias_hit <- intersect(scenario, names(aliases))
  if (length(alias_hit)) {
    warning(
      "Legacy medium aliases were mapped to human presets: ",
      paste(paste0(alias_hit, " -> ", aliases[alias_hit]), collapse = ", "),
      call. = FALSE
    )
    scenario <- unique(c(setdiff(scenario, alias_hit), unname(aliases[alias_hit])))
  }
  if (!is.null(custom_medium) && !is.null(custom_metabolites)) {
    stop("Supply only one of `custom_medium` or `custom_metabolites`.",
         call. = FALSE)
  }

  pieces <- list()
  preset_ids <- intersect(
    scenario,
    c(
      "normal_human_plasma", "high_glucose", "low_glucose",
      "high_lactate", "low_lactate", "rpmi1640"
    )
  )
  for (preset_id in preset_ids) {
    pieces[[length(pieces) + 1L]] <- .rc_build_human_medium_preset(
      gem = gem,
      preset_id = preset_id,
      exchange_limit = exchange_limit,
      uptake_scale = uptake_scale,
      condition = condition %||% "all",
      exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching
    )
  }

  if ("custom" %in% scenario && !is.null(custom_metabolites)) {
    required <- c("metabolite_name", "metabolite_pattern", "available")
    missing <- setdiff(required, colnames(custom_metabolites))
    if (length(missing)) {
      stop("`custom_metabolites` missing columns: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    compounds <- custom_metabolites[
      custom_metabolites$available %in% TRUE, , drop = FALSE
    ]
    if (!nrow(compounds)) {
      stop("`custom_metabolites` contains no available metabolites.",
           call. = FALSE)
    }
    defaults <- list(
      concentration_mM = NA_real_, uptake_fraction = 1,
      target_exchange_flag = FALSE, required_match = TRUE
    )
    for (name in names(defaults)) {
      if (!name %in% colnames(compounds)) compounds[[name]] <- defaults[[name]]
    }
    custom_reference <- data.frame(
      preset_id = "custom",
      reference_label = if ("reference_label" %in% colnames(compounds)) {
        paste(unique(na.omit(as.character(compounds$reference_label))),
              collapse = ";")
      } else {
        "user supplied"
      },
      reference_doi = if ("reference_doi" %in% colnames(compounds)) {
        paste(unique(na.omit(as.character(compounds$reference_doi))),
              collapse = ";")
      } else {
        NA_character_
      },
      reference_pmid = if ("reference_pmid" %in% colnames(compounds)) {
        paste(unique(na.omit(as.character(compounds$reference_pmid))),
              collapse = ";")
      } else {
        NA_character_
      },
      species = "Homo sapiens",
      evidence_scope = "user-supplied metabolite availability",
      stringsAsFactors = FALSE
    )
    pieces[[length(pieces) + 1L]] <- .rc_build_human_medium_preset(
      gem = gem,
      preset_id = "custom",
      exchange_limit = exchange_limit,
      uptake_scale = uptake_scale,
      condition = condition %||% "all",
      exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching,
      compounds = compounds,
      custom_reference = custom_reference
    )
  }

  legacy <- setdiff(scenario, preset_ids)
  if ("custom" %in% legacy && !is.null(custom_metabolites)) {
    legacy <- setdiff(legacy, "custom")
  }
  if (length(legacy)) {
    pieces[[length(pieces) + 1L]] <-
      .rc_medium_presets_previous_make_medium_scenarios(
        gem = gem,
        scenario = legacy,
        custom_medium = custom_medium,
        uptake_scale = uptake_scale,
        condition_col = condition_col,
        exchange_roles = exchange_roles,
        condition = condition,
        exchange_limit = exchange_limit
      )
  }

  output <- .rc_bind_frames_fill(pieces)
  if (!nrow(output)) stop("No medium rows were produced.", call. = FALSE)
  output
}
