# Literature-backed extracellular-medium contract.
#
# R/medium.R retains the low-level exchange mapping and custom-medium machinery.
# This file defines the public biological scenarios and their evidence policy.
# Nutrient composition is taken from high-authority, directly reproducible
# formulations or quantitative extracellular metabolomics. Challenge papers
# supply only the named treatment concentration; they do not define the basal
# nutrient composition unless a complete formulation is provided.

.rc_make_medium_scenarios_unrestricted <- rc_make_medium_scenarios

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

  # Keep only components with an exact HPLM concentration. Rounded serum-ion
  # substitutions and availability-only additions are excluded rather than
  # merged into the authoritative culture background.
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

  # Conservative subset explicitly supported by the Nature 2026 plasma and
  # tissue-fluid metabolomics study. The study quantified 124 metabolites in
  # plasma and extracellular fluids across NSG and C57BL/6J mice. Components
  # outside this auditable set are omitted rather than inherited from HPLM.
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
      target = "glucose", concentration_mM = 25, reference_high_mM = 25,
      concentration_basis = "Han_2015_high_glucose_25mM",
      challenge_reference_label =
        "Han et al., Gynecologic Oncology 2015; 25 mM glucose",
      challenge_reference_doi = "10.1016/j.ygyno.2015.06.036"
    ),
    low_glucose = list(
      target = "glucose", concentration_mM = 1, reference_high_mM = 25,
      concentration_basis = "Han_2015_low_glucose_1mM",
      challenge_reference_label =
        "Han et al., Gynecologic Oncology 2015; 1 mM glucose",
      challenge_reference_doi = "10.1016/j.ygyno.2015.06.036"
    ),
    high_lactate = list(
      target = "lactate", concentration_mM = 20, reference_high_mM = 20,
      concentration_basis = "San_Millan_2020_high_lactate_20mM",
      challenge_reference_label =
        "San-Millan et al., Frontiers in Oncology 2020; 20 mM lactate",
      challenge_reference_doi = "10.3389/fonc.2019.01536"
    ),
    low_lactate = list(
      target = "lactate", concentration_mM = 0.5, reference_high_mM = 20,
      concentration_basis = "Cho_2025_low_lactate_0.5mM",
      challenge_reference_label =
        "Cho et al., Physiological Reports 2025; 0.5 mM lactate",
      challenge_reference_doi = "10.14814/phy2.70450"
    ),
    low_glutamine = list(
      target = "glutamine", concentration_mM = 0.5, reference_high_mM = 4,
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
  compounds$uptake_fraction[selected] <-
    definition$concentration_mM / definition$reference_high_mM
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
      "the named nutrient concentration overridden by the challenge paper."
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
  out$medium_background_id <- background$background_id
  out$background_reference_label <- background$background_reference_label
  out$background_reference_doi <- background$background_reference_doi
  out$background_validation_reference_label <-
    background$background_validation_reference_label
  out$background_validation_reference_doi <-
    background$background_validation_reference_doi
  out$challenge_reference_label <- definition$challenge_reference_label
  out$challenge_reference_doi <- definition$challenge_reference_doi
  out$scenario_construction <-
    "authoritative_HPLM_background_plus_named_nutrient_override"
  out
}

#' Build literature-backed extracellular medium scenarios
#'
#' Human plasma and culture backgrounds use exact HPLM composition from Cell
#' 2017 and the updated formulation in Cell Metabolism 2021. Plasmax from Science
#' Advances 2019 is an independent validation source and is not numerically
#' averaged with HPLM. Mouse plasma uses a conservative availability set anchored
#' to absolute metabolite measurements in Nature 2026; unsupported components are
#' omitted. Challenge papers provide only the named nutrient concentration.
#'
#' @param gem A validated RegCompass GEM.
#' @param scenario One or more identifiers from `"normal_human_plasma"`,
#'   `"mouse_plasma"`, `"high_glucose"`, `"low_glucose"`,
#'   `"high_lactate"`, `"low_lactate"`, `"low_glutamine"`, and
#'   `"custom"`. `NULL` is also accepted for a custom-only run.
#' @param species `"auto"`, `"human"`, or `"mouse"`. Human plasma and culture
#'   challenges require Human-GEM; mouse plasma requires Mouse-GEM.
#' @param custom_medium User-defined reaction-level bounds. Publication columns
#'   are optional but retained when supplied.
#' @param custom_metabolites User-defined metabolite availability and optional
#'   concentration rows. Publication columns are optional but retained.
#' @param uptake_scale Non-negative global or named sensitivity multipliers for
#'   target nutrient relative caps. These are modelling assumptions, not measured
#'   transporter rates.
#' @param exchange_roles Reaction roles treated as exchange reactions.
#' @param condition Shared condition label; canonical scoring requires `"all"`.
#' @param exchange_limit Shared modelling cap intersected with original GEM
#'   directionality.
#' @param strict_preset_matching Stop when required components cannot be mapped
#'   one-to-one to GEM exchanges.
#' @return A reaction-level medium table with composition and challenge citations.
#' @export
rc_make_medium_scenarios <- function(
    gem,
    scenario = "normal_human_plasma",
    species = c("auto", "human", "mouse"),
    custom_medium = NULL,
    custom_metabolites = NULL,
    uptake_scale = 1,
    exchange_roles = c("exchange"),
    condition = "all",
    exchange_limit = 1,
    strict_preset_matching = TRUE) {
  species <- .rc_infer_gem_species(gem, species)
  if (!is.null(custom_medium) && !is.null(custom_metabolites)) {
    stop("Supply only one of `custom_medium` or `custom_metabolites`.",
         call. = FALSE)
  }
  custom_supplied <- !is.null(custom_medium) || !is.null(custom_metabolites)
  if (is.null(scenario)) {
    scenario <- character()
  } else {
    scenario <- unique(trimws(as.character(scenario)))
    if (anyNA(scenario) || any(!nzchar(scenario))) {
      stop("`scenario` must contain non-empty identifiers or be NULL.",
           call. = FALSE)
    }
  }

  choices <- c(
    "normal_human_plasma", "mouse_plasma",
    "high_glucose", "low_glucose", "high_lactate", "low_lactate",
    "low_glutamine", "custom"
  )
  invalid <- setdiff(scenario, choices)
  if (length(invalid)) {
    stop(
      "Unsupported biological medium scenario: ",
      paste(invalid, collapse = ", "),
      ". Technical GEM boundary modes and incomplete presets are not accepted.",
      call. = FALSE
    )
  }
  custom_requested <- "custom" %in% scenario || custom_supplied
  scenario <- setdiff(scenario, "custom")
  if (custom_requested && !custom_supplied) {
    stop(
      "`custom_medium` or `custom_metabolites` is required for `scenario = 'custom'`.",
      call. = FALSE
    )
  }
  if (!length(scenario) && !custom_supplied) {
    stop("At least one built-in or user-defined medium is required.",
         call. = FALSE)
  }

  human_only <- c(
    "normal_human_plasma", "high_glucose", "low_glucose",
    "high_lactate", "low_lactate", "low_glutamine"
  )
  invalid_human <- intersect(scenario, human_only)
  if (identical(species, "mouse") && length(invalid_human)) {
    stop(
      "Human-derived medium scenarios cannot be used with Mouse-GEM: ",
      paste(invalid_human, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (identical(species, "human") && "mouse_plasma" %in% scenario) {
    stop("`mouse_plasma` requires Mouse-GEM.", call. = FALSE)
  }

  pieces <- list()
  if ("normal_human_plasma" %in% scenario) {
    pieces[[length(pieces) + 1L]] <- .rc_build_authoritative_human_plasma(
      gem = gem,
      condition = condition,
      exchange_limit = exchange_limit,
      uptake_scale = uptake_scale,
      exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching
    )
  }
  if ("mouse_plasma" %in% scenario) {
    pieces[[length(pieces) + 1L]] <- .rc_build_authoritative_mouse_plasma(
      gem = gem,
      condition = condition,
      exchange_limit = exchange_limit,
      uptake_scale = uptake_scale,
      exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching
    )
  }

  challenge_ids <- intersect(
    scenario,
    c(
      "high_glucose", "low_glucose", "high_lactate", "low_lactate",
      "low_glutamine"
    )
  )
  for (scenario_id in challenge_ids) {
    pieces[[length(pieces) + 1L]] <- .rc_build_literature_challenge(
      gem = gem,
      scenario_id = scenario_id,
      condition = condition,
      exchange_limit = exchange_limit,
      uptake_scale = uptake_scale,
      exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching
    )
  }

  if (!is.null(custom_medium)) {
    pieces[[length(pieces) + 1L]] <-
      .rc_make_medium_scenarios_unrestricted(
        gem = gem,
        scenario = "custom",
        species = species,
        custom_medium = custom_medium,
        uptake_scale = uptake_scale,
        exchange_roles = exchange_roles,
        condition = condition,
        exchange_limit = exchange_limit,
        strict_preset_matching = strict_preset_matching
      )
  }
  if (!is.null(custom_metabolites)) {
    pieces[[length(pieces) + 1L]] <-
      .rc_make_medium_scenarios_unrestricted(
        gem = gem,
        scenario = "custom",
        species = species,
        custom_metabolites = custom_metabolites,
        uptake_scale = uptake_scale,
        exchange_roles = exchange_roles,
        condition = condition,
        exchange_limit = exchange_limit,
        strict_preset_matching = strict_preset_matching
      )
  }

  diagnostics <- .rc_bind_frames_fill(lapply(
    pieces,
    function(piece) attr(piece, "preset_diagnostics") %||% data.frame()
  ))
  output <- .rc_bind_frames_fill(pieces)
  if (!nrow(output)) stop("No medium rows were produced.", call. = FALSE)
  attr(output, "preset_diagnostics") <- diagnostics
  attr(output, "species") <- species
  attr(output, "medium_policy") <-
    "authoritative_journal_composition_with_explicit_overrides"
  output
}
