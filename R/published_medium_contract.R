# Publication-bound extracellular-medium contract.
#
# The original implementation remains an internal builder for exact exchange
# mapping. The public entry point below restricts built-in scenario identifiers
# to formulations that are explicitly tied to a published paper and prevents
# unpublished technical or synthetic media from being presented as presets.

.rc_make_medium_scenarios_unrestricted <- rc_make_medium_scenarios

.rc_valid_doi <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & grepl("^10\\.[0-9]{4,9}/[^[:space:]]+$", x)
}

.rc_require_published_custom_medium <- function(x, argument) {
  if (!is.data.frame(x) || !nrow(x)) {
    stop("`", argument, "` must be a non-empty data frame.", call. = FALSE)
  }
  required <- c("reference_label", "reference_doi")
  missing <- setdiff(required, colnames(x))
  if (length(missing)) {
    stop(
      "`", argument, "` must include publication provenance columns: ",
      paste(required, collapse = ", "), ".",
      call. = FALSE
    )
  }
  label <- trimws(as.character(x$reference_label))
  doi <- trimws(as.character(x$reference_doi))
  if (anyNA(label) || any(!nzchar(label)) || any(!.rc_valid_doi(doi))) {
    stop(
      "Every custom medium row must have a non-empty `reference_label` and ",
      "a valid published-paper `reference_doi`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_build_cantor2017_hplm <- function(
    gem, condition, exchange_limit, exchange_roles,
    strict_preset_matching) {
  compounds <- .rc_medium_catalog("normal_human_plasma", "human")
  doi <- "10.1016/j.cell.2017.03.023"
  keep <- !is.na(compounds$component_reference_doi) &
    compounds$component_reference_doi == doi &
    is.finite(compounds$concentration_mM)
  compounds <- compounds[keep, , drop = FALSE]
  if (!nrow(compounds)) {
    stop("The Cantor 2017 HPLM component table is unavailable.", call. = FALSE)
  }

  # Published concentrations are retained as provenance. RegCompass does not
  # convert them into transporter rates. Every represented component is treated
  # as available and bounded only by the shared modelling exchange cap and the
  # original GEM directionality.
  compounds$uptake_fraction <- 1
  compounds$target_exchange_flag <- FALSE
  compounds$required_match <- TRUE
  compounds$concentration_basis <- "Cantor_2017_HPLM_formulation"
  compounds$component_reference_doi <- doi

  reference <- data.frame(
    preset_id = "cantor2017_hplm",
    species = "Homo sapiens",
    reference_label = paste(
      "Cantor et al., Cell 2017:",
      "Physiologic Medium Rewires Cellular Metabolism and Reveals Uric Acid",
      "as an Endogenous Inhibitor of UMP Synthase"
    ),
    reference_doi = doi,
    reference_pmid = "28388410",
    evidence_scope = paste(
      "HPLM components with exact published concentrations and direct",
      "one-to-one GEM exchange mapping; ambiguous salt-to-free-ion",
      "conversions are omitted rather than approximated"
    ),
    stringsAsFactors = FALSE
  )

  out <- .rc_build_medium_preset(
    gem = gem,
    preset_id = "cantor2017_hplm",
    species = "human",
    exchange_limit = exchange_limit,
    uptake_scale = 1,
    condition = condition %||% "all",
    exchange_roles = exchange_roles,
    strict_preset_matching = strict_preset_matching,
    compounds = compounds,
    custom_reference = reference
  )
  out$medium_scenario_id <- "cantor2017_hplm"
  out$evidence_source <- "published_HPLM_Cantor_2017"
  out$assumption_level <- "published_component_availability_no_flux_inference"
  out$concentration_used_for_rate_bound <- FALSE
  out$rate_bound_source <-
    "published_availability_intersected_with_original_gem_directionality"
  attr(out, "medium_policy") <- "published_paper_bound_presets_only"
  diagnostics <- attr(out, "preset_diagnostics")
  if (is.data.frame(diagnostics) && nrow(diagnostics)) {
    diagnostics$medium_scenario_id <- "cantor2017_hplm"
    attr(out, "preset_diagnostics") <- diagnostics
  }
  out
}

#' Build publication-bound extracellular medium scenarios
#'
#' Built-in scenario identifiers are restricted to formulations whose encoded
#' components and concentrations can be traced directly to a published paper.
#' The only current built-in preset is `"cantor2017_hplm"` (Cantor et al., Cell
#' 2017; doi:10.1016/j.cell.2017.03.023). Ambiguous or synthetic presets,
#' single-nutrient challenges, manufacturer-only formulations, and technical GEM
#' boundary modes are deliberately not accepted as biological scenarios.
#'
#' Custom reaction-level or metabolite-level environments remain supported, but
#' they are not named presets: set `scenario = NULL` and provide publication
#' provenance through `reference_label` and `reference_doi` columns.
#'
#' @param gem A validated RegCompass GEM.
#' @param scenario Zero or more publication-bound built-in scenario identifiers.
#'   The accepted built-in value is `"cantor2017_hplm"`. Use `NULL` for a
#'   publication-backed custom environment.
#' @param species `"auto"`, `"human"`, or `"mouse"`. The built-in HPLM preset
#'   requires Human-GEM. Mouse analyses must currently provide a published custom
#'   medium because no partial mouse-plasma preset is retained.
#' @param custom_medium Exact reaction-level rows. Every row must include
#'   `reference_label` and a valid published-paper `reference_doi`.
#' @param custom_metabolites Metabolite-availability rows. Every row must include
#'   `reference_label` and a valid published-paper `reference_doi`.
#' @param uptake_scale Retained for API compatibility. Publication-bound runs
#'   require the exact value `1`; arbitrary rescaling would no longer represent
#'   the cited medium.
#' @param exchange_roles Reaction roles treated as exchange reactions.
#' @param condition Shared condition label; canonical scoring requires `"all"`.
#' @param exchange_limit Shared modelling cap intersected with original GEM
#'   directionality. Published concentrations are not converted into uptake flux.
#' @param strict_preset_matching Stop when a published component cannot be mapped
#'   one-to-one to a GEM exchange.
#' @return A reaction-level medium table with publication provenance.
#' @export
rc_make_medium_scenarios <- function(
    gem,
    scenario = "cantor2017_hplm",
    species = c("auto", "human", "mouse"),
    custom_medium = NULL,
    custom_metabolites = NULL,
    uptake_scale = 1,
    exchange_roles = c("exchange"),
    condition = "all",
    exchange_limit = 1,
    strict_preset_matching = TRUE) {
  species <- .rc_infer_gem_species(gem, species)
  if (!is.numeric(uptake_scale) || length(uptake_scale) != 1L ||
      !is.finite(uptake_scale) || uptake_scale != 1) {
    stop(
      "Publication-bound medium scenarios require `uptake_scale = 1`.",
      call. = FALSE
    )
  }
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

  choices <- "cantor2017_hplm"
  invalid <- setdiff(scenario, choices)
  if (length(invalid)) {
    stop(
      "Unsupported or insufficiently documented medium scenario: ",
      paste(invalid, collapse = ", "),
      ". Accepted built-in scenario: `cantor2017_hplm`; use `scenario = NULL` ",
      "with a DOI-cited custom medium for other published formulations.",
      call. = FALSE
    )
  }
  if (!length(scenario) && !custom_supplied) {
    stop(
      "Provide `scenario = 'cantor2017_hplm'` or a DOI-cited custom medium.",
      call. = FALSE
    )
  }
  if (length(scenario) && !identical(species, "human")) {
    stop("`cantor2017_hplm` requires a Human-GEM run.", call. = FALSE)
  }

  pieces <- list()
  if ("cantor2017_hplm" %in% scenario) {
    pieces[[length(pieces) + 1L]] <- .rc_build_cantor2017_hplm(
      gem = gem,
      condition = condition,
      exchange_limit = exchange_limit,
      exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching
    )
  }

  if (!is.null(custom_medium)) {
    .rc_require_published_custom_medium(custom_medium, "custom_medium")
    pieces[[length(pieces) + 1L]] <-
      .rc_make_medium_scenarios_unrestricted(
        gem = gem,
        scenario = "custom",
        species = species,
        custom_medium = custom_medium,
        uptake_scale = 1,
        exchange_roles = exchange_roles,
        condition = condition,
        exchange_limit = exchange_limit,
        strict_preset_matching = strict_preset_matching
      )
  }
  if (!is.null(custom_metabolites)) {
    .rc_require_published_custom_medium(
      custom_metabolites, "custom_metabolites"
    )
    pieces[[length(pieces) + 1L]] <-
      .rc_make_medium_scenarios_unrestricted(
        gem = gem,
        scenario = "custom",
        species = species,
        custom_metabolites = custom_metabolites,
        uptake_scale = 1,
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
  attr(output, "medium_policy") <- "published_paper_bound_presets_only"
  output
}