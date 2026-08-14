# Literature-backed extracellular-medium contract.
#
# R/medium.R retains the low-level exchange mapping and custom-medium machinery.
# This file defines the public biological scenarios and their evidence policy.
# Nutrient composition is taken from high-authority, directly reproducible
# formulations or quantitative extracellular metabolomics. Challenge papers
# supply only the named treatment concentration; they do not define the basal
# nutrient composition unless a complete formulation is provided. Published
# concentrations are provenance and availability evidence, not transporter
# kinetics; they are never converted linearly into exchange flux bounds.

.rc_make_medium_scenarios_unrestricted <- function(
    gem,
    scenario = "custom",
    species = c("auto", "human", "mouse"),
    custom_medium = NULL,
    custom_metabolites = NULL,
    uptake_scale = 1,
    exchange_roles = c("exchange"),
    condition = "all",
    exchange_limit = 1,
    strict_preset_matching = TRUE) {
  species <- .rc_infer_gem_species(gem, species)
  if (!identical(as.character(scenario), "custom")) {
    stop("Internal custom-medium builder accepts only `scenario = 'custom'`.",
         call. = FALSE)
  }
  if (!is.null(custom_medium) && !is.null(custom_metabolites)) {
    stop("Supply only one of `custom_medium` or `custom_metabolites`.",
         call. = FALSE)
  }
  if (is.null(custom_medium) && is.null(custom_metabolites)) {
    stop("A custom reaction-level or metabolite-level medium is required.",
         call. = FALSE)
  }

  if (!is.null(custom_medium)) {
    out <- as.data.frame(custom_medium, stringsAsFactors = FALSE)
    required <- c("exchange_reaction_id", "lb", "ub")
    if (!all(required %in% colnames(out)) || !nrow(out)) {
      stop(
        "`custom_medium` requires non-empty exchange_reaction_id, lb and ub columns.",
        call. = FALSE
      )
    }
    out$exchange_reaction_id <- trimws(as.character(out$exchange_reaction_id))
    out$lb <- suppressWarnings(as.numeric(out$lb))
    out$ub <- suppressWarnings(as.numeric(out$ub))
    if (anyNA(out$exchange_reaction_id) ||
        any(!nzchar(out$exchange_reaction_id)) ||
        anyDuplicated(out$exchange_reaction_id) ||
        any(!is.finite(out$lb)) || any(!is.finite(out$ub))) {
      stop("Custom reaction-level medium rows must be unique and finite.",
           call. = FALSE)
    }
    reactions <- colnames(gem$S)
    if (is.null(reactions) || !length(reactions)) reactions <- names(gem$lb)
    index <- match(out$exchange_reaction_id, reactions)
    if (anyNA(index)) {
      stop(
        "Custom medium contains reaction IDs absent from the GEM: ",
        paste(out$exchange_reaction_id[is.na(index)], collapse = ", "),
        ".", call. = FALSE
      )
    }
    base_lb <- as.numeric(gem$lb)[index]
    base_ub <- as.numeric(gem$ub)[index]
    requested_lb <- out$lb
    requested_ub <- out$ub
    out$lb <- pmax(requested_lb, base_lb)
    out$ub <- pmin(requested_ub, base_ub)
    if (any(out$lb > out$ub)) {
      bad <- out$exchange_reaction_id[out$lb > out$ub]
      stop(
        "Custom medium is incompatible with original GEM directionality for: ",
        paste(bad, collapse = ", "), ".", call. = FALSE
      )
    }
    if (!"available" %in% colnames(out)) out$available <- TRUE
    out$available <- as.logical(out$available)
    if (anyNA(out$available)) {
      stop("Custom medium `available` values must be TRUE or FALSE.",
           call. = FALSE)
    }
    out$lb[!out$available] <- pmax(out$lb[!out$available], 0)
    if (any(out$lb > out$ub)) {
      stop("Unavailable custom exchange rows conflict with base GEM bounds.",
           call. = FALSE)
    }
    if (!"medium_scenario_id" %in% colnames(out)) {
      out$medium_scenario_id <- "custom"
    }
    if (!"condition" %in% colnames(out)) out$condition <- condition %||% "all"
    out$species <- if (identical(species, "human")) {
      "Homo sapiens"
    } else {
      "Mus musculus"
    }
    out$requested_lb <- requested_lb
    out$requested_ub <- requested_ub
    out$original_gem_lb <- base_lb
    out$original_gem_ub <- base_ub
    out$evidence_source <- "user_supplied_custom_medium"
    out$assumption_level <- "explicit_user_flux_bound"
    out$concentration_used_for_rate_bound <- FALSE
    attr(out, "preset_diagnostics") <- data.frame(
      preset_id = unique(as.character(out$medium_scenario_id))[[1L]],
      species = unique(out$species)[[1L]],
      n_requested = nrow(out),
      n_mapped = nrow(out),
      stringsAsFactors = FALSE
    )
    return(out)
  }

  compounds <- as.data.frame(custom_metabolites, stringsAsFactors = FALSE)
  if (!nrow(compounds) || !"metabolite_name" %in% colnames(compounds)) {
    stop("`custom_metabolites` requires a non-empty metabolite_name column.",
         call. = FALSE)
  }
  compounds$metabolite_name <- trimws(as.character(compounds$metabolite_name))
  if (anyNA(compounds$metabolite_name) ||
      any(!nzchar(compounds$metabolite_name)) ||
      anyDuplicated(compounds$metabolite_name)) {
    stop("Custom metabolite names must be unique and non-empty.",
         call. = FALSE)
  }
  if (!"available" %in% colnames(compounds)) compounds$available <- TRUE
  compounds$available <- as.logical(compounds$available)
  if (anyNA(compounds$available)) {
    stop("Custom metabolite `available` values must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (!"concentration_mM" %in% colnames(compounds)) {
    compounds$concentration_mM <- NA_real_
  }
  compounds$concentration_mM <- suppressWarnings(
    as.numeric(compounds$concentration_mM)
  )
  if (!"uptake_fraction" %in% colnames(compounds)) {
    compounds$uptake_fraction <- 1
  }
  compounds$uptake_fraction <- suppressWarnings(
    as.numeric(compounds$uptake_fraction)
  )
  if (any(!is.finite(compounds$uptake_fraction)) ||
      any(compounds$uptake_fraction < 0)) {
    stop("Custom metabolite uptake_fraction must be finite and non-negative.",
         call. = FALSE)
  }
  compounds$uptake_fraction[!compounds$available] <- 0
  if (!"category" %in% colnames(compounds)) compounds$category <- "custom"
  if (!"metabolite_pattern" %in% colnames(compounds)) {
    compounds$metabolite_pattern <- vapply(
      compounds$metabolite_name, .rc_medium_pattern, character(1)
    )
  }
  if (!"gem_metabolite_aliases" %in% colnames(compounds)) {
    compounds$gem_metabolite_aliases <- vapply(
      compounds$metabolite_name,
      function(x) paste(.rc_medium_gem_aliases(x), collapse = ";"),
      character(1)
    )
  }
  if (!"target_exchange_flag" %in% colnames(compounds)) {
    compounds$target_exchange_flag <- TRUE
  }
  if (!"required_match" %in% colnames(compounds)) {
    compounds$required_match <- TRUE
  }
  if (!"concentration_basis" %in% colnames(compounds)) {
    compounds$concentration_basis <- "user_supplied"
  }
  if (!"component_reference_doi" %in% colnames(compounds)) {
    compounds$component_reference_doi <- NA_character_
  }
  reference <- data.frame(
    preset_id = "custom",
    species = if (identical(species, "human")) "Homo sapiens" else "Mus musculus",
    reference_label = "user supplied custom metabolite medium",
    reference_doi = NA_character_,
    reference_pmid = NA_character_,
    evidence_scope = paste(
      "User-supplied metabolite availability and explicit uptake fractions;",
      "concentration values are provenance unless the user separately supplies",
      "a flux sensitivity assumption."
    ),
    stringsAsFactors = FALSE
  )
  out <- .rc_build_medium_preset(
    gem = gem,
    preset_id = "custom",
    species = species,
    exchange_limit = exchange_limit,
    uptake_scale = uptake_scale,
    condition = condition %||% "all",
    exchange_roles = exchange_roles,
    strict_preset_matching = strict_preset_matching,
    compounds = compounds,
    custom_reference = reference
  )
  out$evidence_source <- "user_supplied_custom_medium"
  out$assumption_level <- "explicit_user_metabolite_availability_or_uptake_fraction"
  out$concentration_used_for_rate_bound <- FALSE
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

.rc_medium_concentration_provenance_only <- function(out) {
  if (!is.data.frame(out) || !nrow(out)) return(out)
  out$concentration_used_for_rate_bound <- FALSE
  out$rate_bound_source <-
    "availability_and_explicit_uptake_scale_only_concentration_is_provenance"
  out$assumption_level <-
    "availability_with_uniform_exchange_cap_concentration_provenance_only"
  out
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
      "and is not numerically averaged with HPLM. Concentrations are provenance",
      "and availability evidence and are not converted to exchange rates."
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
  out <- .rc_medium_concentration_provenance_only(out)
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
  compounds$component_reference_doi[quantitative] <- quantitative_secondary
  compounds$uptake_fraction <- 1
  compounds$target_exchange_flag <- FALSE
  compounds$required_match[!quantitative] <- FALSE
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
      "glucose, lactate and glutamine retain quantitative concentration",
      "provenance, but those values are not converted to exchange rates."
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
  out <- .rc_medium_concentration_provenance_only(out)
  out$medium_scenario_id <- "mouse_plasma"
  out$medium_background_id <- "Abbott_2026_Nature_mouse_plasma"
  out$composition_primary_reference_doi <- primary_doi
  out$quantitative_secondary_reference_doi <- quantitative_secondary
  out$scenario_construction <-
    "Nature_2026_mouse_plasma_availability_with_concentration_provenance_only"
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
  # Published challenge concentration is descriptive evidence, not a flux-rate
  # model. Keep the default relative uptake factor at one; an explicit user
  # `uptake_scale` remains the only relative sensitivity assumption here.
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
      "the named nutrient concentration overridden by the challenge paper;",
      "the concentration is provenance and is not converted to a flux bound."
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
  out <- .rc_medium_concentration_provenance_only(out)
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
    "authoritative_HPLM_background_plus_named_concentration_provenance"
  out
}
