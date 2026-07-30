# Literature-backed extracellular-medium contract.
#
# R/medium.R retains the low-level mapping and custom-medium implementation. The
# public entry point below exposes only biological scenarios supported by
# published extracellular measurements or published cell-culture experiments.
# Technical GEM boundary modes remain internal and are not biological media.

.rc_make_medium_scenarios_unrestricted <- rc_make_medium_scenarios

.rc_valid_doi <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & grepl("^10\\.[0-9]{4,9}/[^[:space:]]+$", x)
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

.rc_culture_union_background <- function() {
  rpm <- .rc_medium_catalog("rpmi1640", "human")
  dmem <- .rc_medium_catalog("dmem_high_glucose", "human")
  rpm$background_source <- "RPMI_1640"
  dmem$background_source <- "DMEM"
  combined <- .rc_bind_frames_fill(list(rpm, dmem))
  groups <- split(seq_len(nrow(combined)), combined$metabolite_name)
  rows <- lapply(groups, function(index) {
    block <- combined[index, , drop = FALSE]
    row <- block[1L, , drop = FALSE]
    row$concentration_mM <- NA_real_
    row$uptake_fraction <- 1
    row$target_exchange_flag <- FALSE
    row$required_match <- any(block$required_match %in% TRUE)
    row$concentration_basis <-
      "published_RPMI_DMEM_component_availability_union"
    row$component_reference_doi <- .rc_collapse_nonempty(
      block$component_reference_doi
    )
    row$background_source <- .rc_collapse_nonempty(block$background_source)
    row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.rc_prepare_challenge_background <- function(scenario_id) {
  if (scenario_id %in% c("high_glucose", "low_glucose")) {
    return(list(
      compounds = .rc_culture_union_background(),
      background_id = "published_RPMI_DMEM_nutrient_union",
      background_reference_label = paste(
        "Moore et al., JAMA 1967 RPMI lineage;",
        "Dulbecco and Freeman, Virology 1959 DMEM lineage"
      ),
      background_reference_doi = paste(
        "10.1001/jama.1967.03120080053007",
        "10.1016/0042-6822(59)90063-3",
        sep = ";"
      )
    ))
  }
  if (identical(scenario_id, "high_lactate")) {
    compounds <- .rc_medium_catalog("dmem_high_glucose", "human")
    compounds$uptake_fraction <- 1
    compounds$target_exchange_flag <- FALSE
    return(list(
      compounds = compounds,
      background_id = "published_DMEM_nutrient_background",
      background_reference_label =
        "Dulbecco and Freeman, Virology 1959; DMEM nutrient formulation",
      background_reference_doi = "10.1016/0042-6822(59)90063-3"
    ))
  }
  if (identical(scenario_id, "low_lactate")) {
    compounds <- .rc_medium_catalog("normal_human_plasma", "human")
    compounds$uptake_fraction <- 1
    compounds$target_exchange_flag <- FALSE
    return(list(
      compounds = compounds,
      background_id = "published_plasma_like_nutrient_background",
      background_reference_label = paste(
        "Cantor et al., Cell 2017 HPLM;",
        "Vande Voorde et al., Science Advances 2019 Plasmax"
      ),
      background_reference_doi = paste(
        "10.1016/j.cell.2017.03.023",
        "10.1126/sciadv.aau7314",
        sep = ";"
      )
    ))
  }
  if (identical(scenario_id, "low_glutamine")) {
    compounds <- .rc_medium_catalog("dmem_high_glucose", "human")
    compounds$uptake_fraction <- 1
    compounds$target_exchange_flag <- FALSE
    return(list(
      compounds = compounds,
      background_id = "published_DMEM_nutrient_background",
      background_reference_label = paste(
        "Dulbecco and Freeman, Virology 1959 DMEM lineage;",
        "Visagie et al., Cell Bioscience 2015 DMEM deprivation system"
      ),
      background_reference_doi = paste(
        "10.1016/0042-6822(59)90063-3",
        "10.1186/s13578-015-0030-1",
        sep = ";"
      )
    ))
  }
  stop("Unsupported culture challenge scenario: ", scenario_id, call. = FALSE)
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
    stop("Unsupported culture challenge scenario: ", scenario_id, call. = FALSE)
  )
}

.rc_build_literature_challenge <- function(
    gem, scenario_id, condition, exchange_limit, uptake_scale,
    exchange_roles, strict_preset_matching) {
  background <- .rc_prepare_challenge_background(scenario_id)
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
    definition$challenge_reference_doi
  ))
  reference <- data.frame(
    preset_id = scenario_id,
    species = "Homo sapiens",
    reference_label = paste(
      background$background_reference_label,
      definition$challenge_reference_label,
      sep = "; "
    ),
    reference_doi = all_doi,
    reference_pmid = NA_character_,
    evidence_scope = paste(
      "Composite literature-backed cell-culture environment:",
      background$background_id,
      "with the named nutrient concentration overridden by the challenge paper."
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
  out$challenge_reference_label <- definition$challenge_reference_label
  out$challenge_reference_doi <- definition$challenge_reference_doi
  out$scenario_construction <-
    "published_background_plus_named_nutrient_override"
  out
}

#' Build literature-backed extracellular medium scenarios
#'
#' Biological presets combine published plasma measurements or published cell
#' culture formulations with explicit challenge concentrations. Plasma scenarios
#' may synthesize several published sources to represent nutrients known to be
#' available; unsupported quantitative values remain availability-only. Culture
#' challenges retain the usual nutrients of their published basal medium and
#' override only the named glucose, lactate, or glutamine concentration.
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
#' @return A reaction-level medium table with background and challenge citations.
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
  if (identical(species, "mouse") && length(intersect(scenario, human_only))) {
    stop(
      "Human-derived medium scenarios cannot be used with Mouse-GEM: ",
      paste(intersect(scenario, human_only), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (identical(species, "human") && "mouse_plasma" %in% scenario) {
    stop("`mouse_plasma` requires Mouse-GEM.", call. = FALSE)
  }

  pieces <- list()
  plasma_ids <- intersect(
    scenario,
    c("normal_human_plasma", "mouse_plasma")
  )
  if (length(plasma_ids)) {
    pieces[[length(pieces) + 1L]] <-
      .rc_make_medium_scenarios_unrestricted(
        gem = gem,
        scenario = plasma_ids,
        species = species,
        uptake_scale = uptake_scale,
        exchange_roles = exchange_roles,
        condition = condition,
        exchange_limit = exchange_limit,
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
    "published_plasma_or_culture_background_with_explicit_overrides"
  output
}