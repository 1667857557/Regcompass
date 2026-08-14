# Literature-backed extracellular-medium contract.
#
# R/medium.R owns the public rc_make_medium_scenarios() entry point and the
# low-level exchange mapping machinery. This file provides biological preset
# builders plus one explicit custom-medium builder. It never saves, aliases, or
# redefines the public function.

.rc_make_medium_scenarios_unrestricted <- function(
    gem, scenario = "custom", species = c("auto", "human", "mouse"),
    custom_medium = NULL, custom_metabolites = NULL, uptake_scale = 1,
    exchange_roles = c("exchange"), condition = "all", exchange_limit = 1,
    strict_preset_matching = TRUE) {
  species <- .rc_infer_gem_species(gem, species)
  if (!identical(as.character(scenario), "custom")) {
    stop("The internal custom-medium builder accepts only `scenario = 'custom'`.",
         call. = FALSE)
  }
  if (!is.null(custom_medium) && !is.null(custom_metabolites)) {
    stop("Supply only one of `custom_medium` or `custom_metabolites`.",
         call. = FALSE)
  }
  if (is.null(custom_medium) && is.null(custom_metabolites)) {
    stop("A custom medium requires `custom_medium` or `custom_metabolites`.",
         call. = FALSE)
  }

  if (!is.null(custom_medium)) {
    if (!is.data.frame(custom_medium)) {
      stop("`custom_medium` must be a data.frame.", call. = FALSE)
    }
    out <- as.data.frame(custom_medium, stringsAsFactors = FALSE)
    required <- c(
      "medium_scenario_id", "exchange_reaction_id", "lb", "ub", "available"
    )
    missing <- setdiff(required, colnames(out))
    if (length(missing)) {
      stop("`custom_medium` missing columns: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    out$medium_scenario_id <- trimws(as.character(out$medium_scenario_id))
    out$exchange_reaction_id <- trimws(as.character(out$exchange_reaction_id))
    out$lb <- suppressWarnings(as.numeric(out$lb))
    out$ub <- suppressWarnings(as.numeric(out$ub))
    out$available <- as.logical(out$available)
    out$condition <- if ("condition" %in% colnames(out)) {
      trimws(as.character(out$condition))
    } else {
      as.character(condition %||% "all")
    }
    out$condition[is.na(out$condition) | !nzchar(out$condition)] <- "all"
    if (anyNA(out$medium_scenario_id) || any(!nzchar(out$medium_scenario_id)) ||
        anyNA(out$exchange_reaction_id) || any(!nzchar(out$exchange_reaction_id)) ||
        any(!is.finite(out$lb)) || any(!is.finite(out$ub)) ||
        any(out$lb > out$ub) || anyNA(out$available)) {
      stop("Custom reaction-level rows require non-empty IDs, logical availability, and finite ordered bounds.",
           call. = FALSE)
    }
    key <- paste(out$medium_scenario_id, out$exchange_reaction_id,
                 out$condition, sep = "\001")
    if (anyDuplicated(key)) {
      stop("`custom_medium` contains duplicated scenario/reaction/condition rows.",
           call. = FALSE)
    }
    validated <- rc_validate_gem(gem)
    unknown <- setdiff(out$exchange_reaction_id, validated$reactions)
    if (length(unknown) && isTRUE(strict_preset_matching)) {
      stop("Custom medium reactions missing from GEM: ",
           paste(utils::head(unknown, 10L), collapse = ", "), call. = FALSE)
    }
    if (length(unknown)) {
      warning("Dropping custom medium reactions missing from GEM: ",
              paste(utils::head(unknown, 10L), collapse = ", "), call. = FALSE)
      out <- out[out$exchange_reaction_id %in% validated$reactions, , drop = FALSE]
    }
    if (!nrow(out)) stop("No custom reaction-level medium rows remain.", call. = FALSE)
    defaults <- list(
      metabolite_id = NA_character_,
      evidence_source = "user_supplied_custom_medium",
      assumption_level = "user_supplied_explicit_reaction_bounds",
      target_exchange_flag = FALSE,
      concentration_used_for_rate_bound = FALSE,
      rate_bound_source = "user_supplied_explicit_reaction_bounds",
      reference_label = NA_character_, reference_doi = NA_character_,
      reference_pmid = NA_character_,
      evidence_scope = "user-supplied extracellular environment"
    )
    for (field in names(defaults)) {
      if (!field %in% colnames(out)) out[[field]] <- defaults[[field]]
    }
    out$species <- if (identical(species, "human")) "Homo sapiens" else "Mus musculus"
    attr(out, "preset_diagnostics") <- data.frame()
    rownames(out) <- NULL
    return(out)
  }

  if (!is.data.frame(custom_metabolites)) {
    stop("`custom_metabolites` must be a data.frame.", call. = FALSE)
  }
  compounds <- as.data.frame(custom_metabolites, stringsAsFactors = FALSE)
  if (!"metabolite_name" %in% colnames(compounds)) {
    stop("`custom_metabolites` requires `metabolite_name`.", call. = FALSE)
  }
  compounds$metabolite_name <- trimws(as.character(compounds$metabolite_name))
  if (anyNA(compounds$metabolite_name) || any(!nzchar(compounds$metabolite_name))) {
    stop("Custom metabolite names must be non-empty.", call. = FALSE)
  }
  explicit_uptake <- "uptake_fraction" %in% colnames(compounds)
  defaults <- list(
    medium_scenario_id = "custom",
    available = TRUE,
    concentration_mM = NA_real_,
    uptake_fraction = 1,
    category = "custom_metabolite",
    target_exchange_flag = FALSE,
    required_match = TRUE,
    concentration_basis = "user_supplied",
    component_reference_doi = NA_character_,
    reference_label = NA_character_,
    reference_doi = NA_character_
  )
  for (field in names(defaults)) {
    if (!field %in% colnames(compounds)) compounds[[field]] <- defaults[[field]]
  }
  compounds$medium_scenario_id <- trimws(as.character(compounds$medium_scenario_id))
  compounds$available <- as.logical(compounds$available)
  compounds$concentration_mM <- suppressWarnings(as.numeric(compounds$concentration_mM))
  compounds$uptake_fraction <- suppressWarnings(as.numeric(compounds$uptake_fraction))
  compounds$target_exchange_flag <- as.logical(compounds$target_exchange_flag)
  compounds$required_match <- as.logical(compounds$required_match)
  if (anyNA(compounds$medium_scenario_id) || any(!nzchar(compounds$medium_scenario_id)) ||
      anyNA(compounds$available) || anyNA(compounds$uptake_fraction) ||
      any(!is.finite(compounds$uptake_fraction)) || any(compounds$uptake_fraction < 0) ||
      anyNA(compounds$target_exchange_flag) || anyNA(compounds$required_match)) {
    stop("Custom metabolite rows contain invalid scenario, availability, uptake, or matching values.",
         call. = FALSE)
  }
  finite_concentration <- !is.na(compounds$concentration_mM)
  if (any(!is.finite(compounds$concentration_mM[finite_concentration])) ||
      any(compounds$concentration_mM[finite_concentration] < 0)) {
    stop("Custom metabolite concentrations must be missing or finite non-negative values.",
         call. = FALSE)
  }
  if (!"metabolite_pattern" %in% colnames(compounds)) {
    compounds$metabolite_pattern <- vapply(
      compounds$metabolite_name, .rc_medium_pattern, character(1)
    )
  }
  if (!"gem_metabolite_aliases" %in% colnames(compounds)) {
    compounds$gem_metabolite_aliases <- vapply(
      compounds$metabolite_name,
      function(name) paste(.rc_medium_gem_aliases(name), collapse = ";"),
      character(1)
    )
  }

  scenario_ids <- unique(compounds$medium_scenario_id)
  pieces <- lapply(scenario_ids, function(id) {
    one <- compounds[compounds$medium_scenario_id == id, , drop = FALSE]
    reference_label <- .rc_collapse_nonempty(one$reference_label)
    reference_doi <- .rc_collapse_nonempty(one$reference_doi)
    reference <- data.frame(
      preset_id = id,
      species = if (identical(species, "human")) "Homo sapiens" else "Mus musculus",
      reference_label = reference_label,
      reference_doi = reference_doi,
      reference_pmid = NA_character_,
      evidence_scope = "user-supplied extracellular metabolite environment",
      stringsAsFactors = FALSE
    )
    built <- .rc_build_medium_preset(
      gem = gem, preset_id = id, species = species,
      exchange_limit = exchange_limit, uptake_scale = uptake_scale,
      condition = condition %||% "all", exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching,
      compounds = one, custom_reference = reference
    )
    source_index <- match(built$preset_metabolite, one$metabolite_name)
    built$available <- one$available[source_index]
    built$evidence_source <- "user_supplied_custom_metabolites"
    built$assumption_level <- if (explicit_uptake) {
      "user_supplied_relative_uptake_assumption"
    } else {
      "user_supplied_availability_without_flux_conversion"
    }
    built$concentration_used_for_rate_bound <- explicit_uptake &
      built$target_exchange_flag %in% TRUE
    built$rate_bound_source <- ifelse(
      built$concentration_used_for_rate_bound,
      "user_supplied_uptake_fraction_not_inferred_from_concentration",
      "user_supplied_availability_intersected_with_original_gem_directionality"
    )
    built
  })
  diagnostics <- .rc_bind_frames_fill(lapply(
    pieces, function(piece) attr(piece, "preset_diagnostics") %||% data.frame()
  ))
  out <- .rc_bind_frames_fill(pieces)
  attr(out, "preset_diagnostics") <- diagnostics
  rownames(out) <- NULL
  out
}

.rc_collapse_nonempty <- function(x) {
  x <- trimws(as.character(x))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x)) paste(x, collapse = ";") else NA_character_
}

.rc_add_compound_if_missing <- function(compounds, metabolite_name) {
  if (metabolite_name %in% compounds$metabolite_name) return(compounds)
  row <- .rc_medium_rows(
    metabolite_name,
    concentration_mM = NA_real_,
    category = "challenge_nutrient",
    required = TRUE
  )
  row$concentration_basis <- "challenge_target_only"
  row$component_reference_doi <- NA_character_
  .rc_bind_frames_fill(list(compounds, row))
}

.rc_authoritative_hplm_background <- function() {
  compounds <- .rc_medium_catalog("normal_human_plasma", "human")
  cell_2017 <- "10.1016/j.cell.2017.03.023"
  cell_metabolism_2021 <- "10.1016/j.cmet.2021.02.005"
  updated_components <- c(
    "alpha_ketoglutarate", "acetylcarnitine", "malate", "uridine"
  )

  keep <- is.finite(compounds$concentration_mM) & (
    compounds$component_reference_doi == cell_2017 |
      compounds$metabolite_name %in% updated_components
  )
  keep[is.na(keep)] <- FALSE
  compounds <- compounds[keep, , drop = FALSE]
  if (!nrow(compounds)) {
    stop("The authoritative HPLM component table is unavailable.", call. = FALSE)
  }

  updated <- compounds$metabolite_name %in% updated_components
  compounds$component_reference_doi[!updated] <- cell_2017
  compounds$component_reference_doi[updated] <- cell_metabolism_2021
  compounds$concentration_basis[!updated] <-
    "Cantor_2017_Cell_HPLM_formulation"
  compounds$concentration_basis[updated] <-
    "Rossiter_2021_Cell_Metabolism_updated_HPLM_formulation"
  compounds$uptake_fraction <- 1
  compounds$target_exchange_flag <- FALSE
  compounds$background_source <- ifelse(
    updated,
    "Cell_Metabolism_2021_updated_HPLM",
    "Cell_2017_HPLM"
  )
  rownames(compounds) <- NULL
  compounds
}

.rc_build_authoritative_human_plasma <- function(
    gem, condition, exchange_limit, uptake_scale,
    exchange_roles, strict_preset_matching) {
  primary_doi <- paste(
    "10.1016/j.cell.2017.03.023",
    "10.1016/j.cmet.2021.02.005",
    sep = ";"
  )
  validation_doi <- "10.1126/sciadv.aau7314"
  reference <- data.frame(
    preset_id = "normal_human_plasma",
    species = "Homo sapiens",
    reference_label = paste(
      "Cantor et al., Cell 2017 HPLM;",
      "Rossiter et al., Cell Metabolism 2021 updated HPLM;",
      "Vande Voorde et al., Science Advances 2019 Plasmax validation"
    ),
    reference_doi = paste(primary_doi, validation_doi, sep = ";"),
    reference_pmid = "28388410;33651980;30613774",
    evidence_scope = paste(
      "Exact HPLM component concentrations from Cell and Cell Metabolism;",
      "Plasmax is retained as an independent physiological-medium validation",
      "and is not numerically averaged with HPLM."
    ),
    stringsAsFactors = FALSE
  )

  out <- .rc_build_medium_preset(
    gem = gem,
    preset_id = "normal_human_plasma",
    species = "human",
    exchange_limit = exchange_limit,
    uptake_scale = uptake_scale,
    condition = condition %||% "all",
    exchange_roles = exchange_roles,
    strict_preset_matching = strict_preset_matching,
    compounds = .rc_authoritative_hplm_background(),
    custom_reference = reference
  )
  out$medium_scenario_id <- "normal_human_plasma"
  out$medium_background_id <- "authoritative_HPLM_2017_2021"
  out$composition_primary_reference_doi <- primary_doi
  out$composition_validation_reference_doi <- validation_doi
  out$scenario_construction <-
    "authoritative_HPLM_composition_without_cross_study_averaging"
  out
}

.rc_authoritative_mouse_compounds <- function() {
  compounds <- .rc_medium_catalog("mouse_plasma", "mouse")
  nature_2026 <- "10.1038/s41586-025-09898-9"
  quantitative_secondary <- "10.1152/ajpcell.00452.2024"
  supported <- c(
    "glucose", "lactate", "glutamine", "arginine", "ornithine",
    "citrulline", "isoleucine", "leucine", "valine", "serine",
    "asparagine", "proline", "hypoxanthine", "uridine", "creatine"
  )
  compounds <- compounds[
    compounds$metabolite_name %in% supported,
    ,
    drop = FALSE
  ]
  if (!nrow(compounds)) {
    stop("The authoritative mouse-plasma component table is unavailable.",
         call. = FALSE)
  }

  quantitative <- compounds$metabolite_name %in%
    c("glucose", "lactate", "glutamine")
  compounds$component_reference_doi[!quantitative] <- nature_2026
  compounds$concentration_basis[!quantitative] <-
    "Abbott_2026_Nature_mouse_plasma_detected_availability"
  compounds$uptake_fraction[!quantitative] <- 1
  compounds$target_exchange_flag[!quantitative] <- FALSE
  compounds$required_match[!quantitative] <- FALSE
  compounds$component_reference_doi[quantitative] <- quantitative_secondary
  rownames(compounds) <- NULL
  compounds
}

.rc_build_authoritative_mouse_plasma <- function(
    gem, condition, exchange_limit, uptake_scale,
    exchange_roles, strict_preset_matching) {
  primary_doi <- "10.1038/s41586-025-09898-9"
  quantitative_secondary <- "10.1152/ajpcell.00452.2024"
  reference <- data.frame(
    preset_id = "mouse_plasma",
    species = "Mus musculus",
    reference_label = paste(
      "Abbott et al., Nature 2026 absolute metabolite quantification in",
      "mouse plasma and tissue interstitial fluids; Gardner and Stuart 2024",
      "for the retained glucose, lactate and glutamine quantitative values"
    ),
    reference_doi = paste(primary_doi, quantitative_secondary, sep = ";"),
    reference_pmid = NA_character_,
    evidence_scope = paste(
      "Conservative mouse-plasma availability set anchored to Nature 2026;",
      "only glucose, lactate and glutamine retain the secondary quantitative",
      "mouse-plasma values. Unsupported components are omitted."
    ),
    stringsAsFactors = FALSE
  )

  out <- .rc_build_medium_preset(
    gem = gem,
    preset_id = "mouse_plasma",
    species = "mouse",
    exchange_limit = exchange_limit,
    uptake_scale = uptake_scale,
    condition = condition %||% "all",
    exchange_roles = exchange_roles,
    strict_preset_matching = strict_preset_matching,
    compounds = .rc_authoritative_mouse_compounds(),
    custom_reference = reference
  )
  out$medium_scenario_id <- "mouse_plasma"
  out$medium_background_id <- "Abbott_2026_Nature_mouse_plasma"
  out$composition_primary_reference_doi <- primary_doi
  out$quantitative_secondary_reference_doi <- quantitative_secondary
  out$scenario_construction <-
    "Nature_2026_mouse_plasma_availability_with_limited_quantitative_secondary"
  out
}

.rc_prepare_challenge_background <- function() {
  list(
    compounds = .rc_authoritative_hplm_background(),
    background_id = "authoritative_HPLM_2017_2021",
    background_reference_label = paste(
      "Cantor et al., Cell 2017 HPLM;",
      "Rossiter et al., Cell Metabolism 2021 updated HPLM"
    ),
    background_reference_doi = paste(
      "10.1016/j.cell.2017.03.023",
      "10.1016/j.cmet.2021.02.005",
      sep = ";"
    ),
    background_validation_reference_label =
      "Vande Voorde et al., Science Advances 2019 Plasmax",
    background_validation_reference_doi = "10.1126/sciadv.aau7314"
  )
}

.rc_challenge_definition <- function(scenario_id) {
  switch(
    scenario_id,
    high_glucose = list(
      target = "glucose", concentration_mM = 25,
      concentration_basis = "Han_2015_high_glucose_25mM",
      challenge_reference_label =
        "Han et al., Gynecologic Oncology 2015; 25 mM glucose",
      challenge_reference_doi = "10.1016/j.ygyno.2015.06.036"
    ),
    low_glucose = list(
      target = "glucose", concentration_mM = 1,
      concentration_basis = "Han_2015_low_glucose_1mM",
      challenge_reference_label =
        "Han et al., Gynecologic Oncology 2015; 1 mM glucose",
      challenge_reference_doi = "10.1016/j.ygyno.2015.06.036"
    ),
    high_lactate = list(
      target = "lactate", concentration_mM = 20,
      concentration_basis = "San_Millan_2020_high_lactate_20mM",
      challenge_reference_label =
        "San-Millan et al., Frontiers in Oncology 2020; 20 mM lactate",
      challenge_reference_doi = "10.3389/fonc.2019.01536"
    ),
    low_lactate = list(
      target = "lactate", concentration_mM = 0.5,
      concentration_basis = "Cho_2025_low_lactate_0.5mM",
      challenge_reference_label =
        "Cho et al., Physiological Reports 2025; 0.5 mM lactate",
      challenge_reference_doi = "10.14814/phy2.70450"
    ),
    low_glutamine = list(
      target = "glutamine", concentration_mM = 0.5,
      concentration_basis = "Visagie_2015_low_glutamine_0.5mM",
      challenge_reference_label = paste(
        "Visagie et al., Cell Bioscience 2015;",
        "Methods-defined 0.5 mM glutamine condition"
      ),
      challenge_reference_doi = "10.1186/s13578-015-0030-1"
    ),
    stop("Unsupported culture challenge scenario: ", scenario_id,
         call. = FALSE)
  )
}

.rc_build_literature_challenge <- function(
    gem, scenario_id, condition, exchange_limit, uptake_scale,
    exchange_roles, strict_preset_matching) {
  background <- .rc_prepare_challenge_background()
  definition <- .rc_challenge_definition(scenario_id)
  compounds <- .rc_add_compound_if_missing(
    background$compounds,
    definition$target
  )
  selected <- compounds$metabolite_name == definition$target
  compounds$concentration_mM[selected] <- definition$concentration_mM
  compounds$concentration_basis[selected] <- definition$concentration_basis
  compounds$component_reference_doi[selected] <-
    definition$challenge_reference_doi
  compounds$uptake_fraction[selected] <- 1
  compounds$target_exchange_flag[selected] <- TRUE
  compounds$required_match[selected] <- TRUE

  all_doi <- .rc_collapse_nonempty(c(
    background$background_reference_doi,
    background$background_validation_reference_doi,
    definition$challenge_reference_doi
  ))
  reference <- data.frame(
    preset_id = scenario_id,
    species = "Homo sapiens",
    reference_label = paste(
      background$background_reference_label,
      background$background_validation_reference_label,
      definition$challenge_reference_label,
      sep = "; "
    ),
    reference_doi = all_doi,
    reference_pmid = NA_character_,
    evidence_scope = paste(
      "Authoritative HPLM basal composition from Cell and Cell Metabolism,",
      "independently validated against Plasmax in Science Advances, with only",
      "the named nutrient concentration overridden as scenario metadata;",
      "concentration is not converted to an uptake flux bound."
    ),
    stringsAsFactors = FALSE
  )

  out <- .rc_build_medium_preset(
    gem = gem,
    preset_id = scenario_id,
    species = "human",
    exchange_limit = exchange_limit,
    uptake_scale = uptake_scale,
    condition = condition %||% "all",
    exchange_roles = exchange_roles,
    strict_preset_matching = strict_preset_matching,
    compounds = compounds,
    custom_reference = reference
  )
  target_rows <- out$target_exchange_flag %in% TRUE
  out$concentration_used_for_rate_bound[target_rows] <- FALSE
  out$rate_bound_source[target_rows] <-
    "background_model_cap_concentration_metadata_only"
  out$assumption_level[target_rows] <-
    "literature_concentration_metadata_without_flux_conversion"
  out$medium_background_id <- background$background_id
  out$background_reference_label <- background$background_reference_label
  out$background_reference_doi <- background$background_reference_doi
  out$background_validation_reference_label <-
    background$background_validation_reference_label
  out$background_validation_reference_doi <-
    background$background_validation_reference_doi
  out$challenge_reference_label <- definition$challenge_reference_label
  out$challenge_reference_doi <- definition$challenge_reference_doi
  out$scenario_construction <- paste(
    "authoritative_HPLM_background_plus_named_nutrient_concentration_metadata;",
    "no_automatic_concentration_to_flux_mapping"
  )
  out
}
