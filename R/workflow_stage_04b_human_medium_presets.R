# Workflow stage 4b: published human medium presets.
#
# Presets control exchange-reaction uptake. Unlisted exchanges remain closed for
# uptake when rc_apply_medium_constraints() applies the returned table; secretion
# remains governed by the existing medium application policy.

.rc_medium_presets_previous_make_medium_scenarios <- rc_make_medium_scenarios

.rc_medium_reference_catalog <- function() {
  data.frame(
    preset_id = c(
      "normal_human_plasma", "high_glucose", "low_glucose",
      "high_lactate", "low_lactate", "rpmi1640"
    ),
    reference_label = c(
      paste(
        "Cantor et al., Cell 2017 (HPLM);",
        "Psychogios et al., PLoS One 2011 (human serum metabolome)"
      ),
      "Han et al., Gynecologic Oncology 2015 (human endometrial cancer cells)",
      "Han et al., Gynecologic Oncology 2015 (human endometrial cancer cells)",
      paste(
        "Schwickert et al., Experientia 1996 (human cervical tumors);",
        "Kennedy et al., PLoS One 2013 (human breast cells/tumors)"
      ),
      "Kennedy et al., PLoS One 2013 (0.5-2 mM human physiologic range)",
      paste(
        "Moore et al., JAMA 1967 (normal human leukocytes);",
        "Cantor et al., Cell 2017 (RPMI versus human plasma)"
      )
    ),
    reference_doi = c(
      "10.1016/j.cell.2017.03.023;10.1371/journal.pone.0016957",
      "10.1016/j.ygyno.2015.06.036",
      "10.1016/j.ygyno.2015.06.036",
      "10.1007/BF01919316;10.1371/journal.pone.0075154",
      "10.1371/journal.pone.0075154",
      "10.1001/jama.1967.03120080053007;10.1016/j.cell.2017.03.023"
    ),
    reference_pmid = c(
      "28388410;21359215", "26135947", "26135947",
      "8641383;24069390", "24069390", "28388410"
    ),
    species = "Homo sapiens",
    evidence_scope = c(
      "adult human plasma/serum polar-metabolite composition",
      "25 mM high-glucose exposure in human cancer cells",
      "1 mM glucose limitation in human cancer cells",
      "20 mM high-lactate human tumor/cell sensitivity condition",
      "0.5 mM lower physiologic human blood lactate condition",
      "serum-free RPMI-1640 basal formulation for human leukocyte culture"
    ),
    stringsAsFactors = FALSE
  )
}

.rc_medium_metabolite_pattern <- function(name) {
  special <- c(
    glucose = "(^|[^a-z0-9])(d[- ]?)?glucose([^a-z0-9]|$)|(^|[^a-z0-9])glc([^a-z0-9]|$)",
    lactate = "(^|[^a-z0-9])(l[- ]?)?lactate([^a-z0-9]|$)|lactic[ _-]?acid",
    pyruvate = "(^|[^a-z0-9])pyruvate([^a-z0-9]|$)|pyruvic[ _-]?acid",
    aspartate = "(^|[^a-z0-9])(l[- ]?)?aspartate([^a-z0-9]|$)|aspartic[ _-]?acid",
    cysteine_cystine = "cysteine|cystine",
    cystine = "cystine|cysteine",
    glutamate = "(^|[^a-z0-9])(l[- ]?)?glutamate([^a-z0-9]|$)|glutamic[ _-]?acid",
    acetate = "(^|[^a-z0-9])acetate([^a-z0-9]|$)|acetic[ _-]?acid",
    citrate = "(^|[^a-z0-9])citrate([^a-z0-9]|$)|citric[ _-]?acid",
    fumarate = "(^|[^a-z0-9])fumarate([^a-z0-9]|$)|fumaric[ _-]?acid",
    malate = "(^|[^a-z0-9])malate([^a-z0-9]|$)|malic[ _-]?acid",
    succinate = "(^|[^a-z0-9])succinate([^a-z0-9]|$)|succinic[ _-]?acid",
    beta_hydroxybutyrate = "beta[- _]?hydroxybutyrate|3[- _]?hydroxybutyrate|b[- _]?hydroxybutyrate",
    acetylcarnitine = "acetyl[- _]?carnitine",
    urate = "(^|[^a-z0-9])urate([^a-z0-9]|$)|uric[ _-]?acid",
    carbon_dioxide = "carbon[ _-]?dioxide|(^|[^a-z0-9])co2([^a-z0-9]|$)",
    bicarbonate = "bicarbonate|(^|[^a-z0-9])hco3([^a-z0-9]|$)",
    phosphate = "orthophosphate|inorganic[ _-]?phosphate|(^|[^a-z0-9])phosphate([^a-z0-9]|$)",
    sulfate = "(^|[^a-z0-9])sulfate([^a-z0-9]|$)|(^|[^a-z0-9])so4([^a-z0-9]|$)",
    oxygen = "(^|[^a-z0-9])oxygen([^a-z0-9]|$)|(^|[^a-z0-9])o2([^a-z0-9]|$)",
    water = "(^|[^a-z0-9])water([^a-z0-9]|$)|(^|[^a-z0-9])h2o([^a-z0-9]|$)",
    sodium = "(^|[^a-z0-9])sodium([^a-z0-9]|$)|(^|[^a-z0-9])na\\+?([^a-z0-9]|$)",
    potassium = "(^|[^a-z0-9])potassium([^a-z0-9]|$)|(^|[^a-z0-9])k\\+?([^a-z0-9]|$)",
    calcium = "(^|[^a-z0-9])calcium([^a-z0-9]|$)|(^|[^a-z0-9])ca2\\+?([^a-z0-9]|$)",
    magnesium = "(^|[^a-z0-9])magnesium([^a-z0-9]|$)|(^|[^a-z0-9])mg2\\+?([^a-z0-9]|$)",
    chloride = "(^|[^a-z0-9])chloride([^a-z0-9]|$)|(^|[^a-z0-9])cl-?([^a-z0-9]|$)",
    ammonium = "(^|[^a-z0-9])ammonium([^a-z0-9]|$)|(^|[^a-z0-9])nh4([^a-z0-9]|$)",
    hydroxyproline = "hydroxyproline",
    inositol = "inositol",
    folate = "folate|folic[ _-]?acid",
    pantothenate = "pantothenate|pantothenic[ _-]?acid",
    niacinamide = "niacinamide|nicotinamide",
    p_aminobenzoate = "p[- _]?aminobenzoate|para[- _]?aminobenzoic",
    pyridoxine = "pyridoxine|vitamin[ _-]?b6",
    riboflavin = "riboflavin|vitamin[ _-]?b2",
    thiamine = "thiamine|vitamin[ _-]?b1",
    cobalamin = "cobalamin|vitamin[ _-]?b12"
  )
  if (name %in% names(special)) return(unname(special[[name]]))
  token <- gsub("_", "[- _]?", name, fixed = TRUE)
  paste0("(^|[^a-z0-9])(l[- ]?)?", token, "([^a-z0-9]|$)")
}

.rc_medium_compound_table <- function(
    names, concentration_mM = NA_real_, uptake_fraction = 1,
    target = FALSE, required = FALSE) {
  n <- length(names)
  recycle <- function(x) rep(x, length.out = n)
  data.frame(
    metabolite_name = names,
    metabolite_pattern = vapply(
      names, .rc_medium_metabolite_pattern, character(1)
    ),
    concentration_mM = recycle(concentration_mM),
    uptake_fraction = recycle(uptake_fraction),
    target_exchange_flag = recycle(target),
    required_match = recycle(required),
    stringsAsFactors = FALSE
  )
}

.rc_medium_background_compounds <- function() {
  .rc_medium_compound_table(
    c(
      "oxygen", "water", "carbon_dioxide", "bicarbonate",
      "sodium", "potassium", "calcium", "magnesium", "chloride",
      "phosphate", "sulfate", "ammonium"
    )
  )
}

.rc_medium_human_plasma_compounds <- function() {
  names <- c(
    "glucose", "lactate", "pyruvate", "alanine", "arginine",
    "asparagine", "aspartate", "cysteine_cystine", "glutamate",
    "glutamine", "glycine", "histidine", "isoleucine", "leucine",
    "lysine", "methionine", "phenylalanine", "proline", "serine",
    "threonine", "tryptophan", "tyrosine", "valine", "acetate",
    "citrate", "fumarate", "malate", "succinate", "glycerol",
    "beta_hydroxybutyrate", "choline", "carnitine",
    "acetylcarnitine", "urate", "uridine"
  )
  concentrations <- c(
    5.0, 1.5, 0.05, 0.43, 0.11,
    0.05, 0.02, 0.04, 0.08,
    0.55, 0.30, 0.11, 0.07, 0.10,
    0.18, 0.03, 0.06, 0.20, 0.12,
    0.14, 0.05, 0.06, 0.20, 0.10,
    0.10, 0.005, 0.005, 0.005, 0.10,
    0.10, 0.01, 0.04, 0.005, 0.35, 0.003
  )
  fractions <- rep(1, length(names))
  fractions[names == "glucose"] <- 5 / 25
  fractions[names == "lactate"] <- 1.5 / 20
  target <- names %in% c("glucose", "lactate")
  required <- names %in% c("glucose", "lactate", "glutamine")
  rbind(
    .rc_medium_compound_table(
      names, concentrations, fractions, target, required
    ),
    .rc_medium_background_compounds()
  )
}

.rc_medium_rpmi1640_compounds <- function() {
  names <- c(
    "glucose", "arginine", "asparagine", "aspartate", "cystine",
    "glutamate", "glutamine", "glycine", "histidine",
    "hydroxyproline", "isoleucine", "leucine", "lysine",
    "methionine", "phenylalanine", "proline", "serine",
    "threonine", "tryptophan", "tyrosine", "valine", "choline",
    "inositol", "folate", "pantothenate", "biotin", "niacinamide",
    "p_aminobenzoate", "pyridoxine", "riboflavin", "thiamine",
    "cobalamin", "glutathione"
  )
  concentrations <- c(
    11.1, 1.149, 0.379, 0.150, 0.208,
    0.136, 2.055, 0.133, 0.097, 0.153,
    0.382, 0.382, 0.219, 0.101, 0.091,
    0.174, 0.286, 0.168, 0.0245, 0.111,
    0.171, 0.0214, 0.194, 0.00227, 0.00052,
    0.00082, 0.0082, 0.0073, 0.0049, 0.00027,
    0.0030, 0.0000037, 0.00325
  )
  fractions <- rep(1, length(names))
  fractions[names == "glucose"] <- 11.1 / 25
  target <- names == "glucose"
  required <- names %in% c("glucose", "glutamine")
  rbind(
    .rc_medium_compound_table(
      names, concentrations, fractions, target, required
    ),
    .rc_medium_background_compounds()
  )
}

.rc_medium_preset_compounds <- function(preset_id) {
  plasma <- .rc_medium_human_plasma_compounds()
  if (identical(preset_id, "normal_human_plasma")) return(plasma)
  if (identical(preset_id, "rpmi1640")) {
    return(.rc_medium_rpmi1640_compounds())
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
  row <- plasma$metabolite_name == target
  plasma$concentration_mM[row] <- concentration
  plasma$uptake_fraction[row] <- concentration / reference_high
  plasma$target_exchange_flag[row] <- TRUE
  plasma$required_match[row] <- TRUE
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

.rc_medium_scenario_scale <- function(uptake_scale, scenario_id) {
  if (!is.numeric(uptake_scale) || !length(uptake_scale) ||
      any(!is.finite(uptake_scale)) || any(uptake_scale < 0)) {
    stop("`uptake_scale` must contain finite non-negative values.",
         call. = FALSE)
  }
  if (length(uptake_scale) == 1L && (
      is.null(names(uptake_scale)) || !nzchar(names(uptake_scale)[[1L]]))) {
    return(as.numeric(uptake_scale[[1L]]))
  }
  if (!is.null(names(uptake_scale)) &&
      scenario_id %in% names(uptake_scale)) {
    return(as.numeric(uptake_scale[[scenario_id]]))
  }
  1
}

.rc_build_human_medium_preset <- function(
    gem, preset_id, exchange_limit, uptake_scale = 1,
    condition = "all", exchange_roles = "exchange",
    strict_preset_matching = TRUE, compounds = NULL,
    custom_reference = NULL) {
  if (!is.numeric(exchange_limit) || length(exchange_limit) != 1L ||
      !is.finite(exchange_limit) || exchange_limit <= 0) {
    stop("`exchange_limit` must be one positive finite number.",
         call. = FALSE)
  }
  if (!is.logical(strict_preset_matching) ||
      length(strict_preset_matching) != 1L ||
      is.na(strict_preset_matching)) {
    stop("`strict_preset_matching` must be TRUE or FALSE.", call. = FALSE)
  }

  validated <- rc_validate_gem(gem)
  if (is.null(gem$reaction_meta) ||
      !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem)
  }
  meta <- gem$reaction_meta[
    match(validated$reactions, as.character(gem$reaction_meta$reaction_id)),
    ,
    drop = FALSE
  ]
  roles <- unique(trimws(as.character(exchange_roles)))
  roles <- roles[!is.na(roles) & nzchar(roles)]
  exchange_meta <- meta[as.character(meta$role) %in% roles, , drop = FALSE]
  exchange_ids <- as.character(exchange_meta$reaction_id)
  if (!length(exchange_ids)) {
    stop("No exchange reactions were identified for the medium preset.",
         call. = FALSE)
  }

  if (is.null(compounds)) compounds <- .rc_medium_preset_compounds(preset_id)
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
  rows <- vector("list", nrow(compounds))
  matched_count <- integer(nrow(compounds))
  scenario_scale <- .rc_medium_scenario_scale(uptake_scale, preset_id)
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
    original_lb <- as.numeric(validated$lb[index])
    original_ub <- as.numeric(validated$ub[index])
    fraction <- min(
      1,
      as.numeric(compounds$uptake_fraction[[i]]) * scenario_scale
    )
    lb <- pmax(original_lb, -exchange_limit * fraction)
    lb[original_lb >= 0] <- 0
    ub <- pmin(original_ub, exchange_limit)
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
      ub = ub,
      available = TRUE,
      original_lb = original_lb,
      original_ub = original_ub,
      exchange_limit = exchange_limit,
      uptake_fraction = fraction,
      evidence_source = "published_human_medium_preset",
      assumption_level =
        "literature_defined_availability_with_relative_uptake_cap",
      target_exchange_flag =
        as.logical(compounds$target_exchange_flag[[i]]),
      concentration_used_for_rate_bound =
        as.logical(compounds$target_exchange_flag[[i]]),
      rate_bound_source = if (
        isTRUE(compounds$target_exchange_flag[[i]])
      ) {
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
  output <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(output)) {
    stop("The medium preset did not match any exchange reactions.",
         call. = FALSE)
  }
  output <- output[
    order(
      !output$target_exchange_flag,
      match(output$preset_metabolite, compounds$metabolite_name)
    ),
    ,
    drop = FALSE
  ]
  duplicate_ids <- unique(
    output$exchange_reaction_id[duplicated(output$exchange_reaction_id)]
  )
  output <- output[
    !duplicated(output$exchange_reaction_id),
    ,
    drop = FALSE
  ]

  reference <- custom_reference %||% .rc_medium_reference_catalog()
  ref <- reference[reference$preset_id == preset_id, , drop = FALSE]
  if (!nrow(ref)) {
    ref <- data.frame(
      preset_id = preset_id,
      reference_label = "user supplied",
      reference_doi = NA_character_,
      reference_pmid = NA_character_,
      species = "Homo sapiens",
      evidence_scope = "user-supplied custom metabolite availability",
      stringsAsFactors = FALSE
    )
  }
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
  attr(output, "preset_references") <- ref
  attr(output, "duplicate_exchange_matches_resolved") <- duplicate_ids
  output
}

rc_make_medium_scenarios <- function(
    gem,
    scenario = "compass_model_bounds",
    custom_medium = NULL,
    custom_metabolites = NULL,
    uptake_scale = 1,
    condition_col = NULL,
    exchange_roles = c("exchange"),
    condition = condition_col,
    exchange_limit = 1,
    strict_preset_matching = TRUE) {
  choices <- c(
    "compass_model_bounds",
    "normal_human_plasma", "high_glucose", "low_glucose",
    "high_lactate", "low_lactate", "rpmi1640",
    "permissive_all_exchange", "minimal", "low_glutamine",
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
      paste(
        paste0(alias_hit, " -> ", unname(aliases[alias_hit])),
        collapse = ", "
      ),
      call. = FALSE
    )
    scenario <- unique(c(
      setdiff(scenario, alias_hit),
      unname(aliases[alias_hit])
    ))
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
      custom_metabolites$available %in% TRUE,
      ,
      drop = FALSE
    ]
    if (!nrow(compounds)) {
      stop("`custom_metabolites` contains no available metabolites.",
           call. = FALSE)
    }
    defaults <- list(
      concentration_mM = NA_real_,
      uptake_fraction = 1,
      target_exchange_flag = FALSE,
      required_match = TRUE
    )
    for (name in names(defaults)) {
      if (!name %in% colnames(compounds)) {
        compounds[[name]] <- defaults[[name]]
      }
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
      evidence_scope = "user-supplied custom metabolite availability",
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
