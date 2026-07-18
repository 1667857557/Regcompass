# Literature-backed extracellular environments --------------------------------

.rc_infer_gem_species <- function(gem, species = c("auto", "human", "mouse")) {
  species <- match.arg(species)
  recorded <- tolower(as.character(gem$model_info$species %||% ""))
  if (identical(species, "auto")) {
    if (recorded %in% c("human", "mouse")) return(recorded)
    source <- tolower(as.character(gem$model_info$source %||% ""))
    if (grepl("mouse", source)) return("mouse")
    if (grepl("human", source)) return("human")
    return("human")
  }
  if (recorded %in% c("human", "mouse") && !identical(recorded, species)) {
    stop(
      "Requested species `", species, "` does not match GEM species `",
      recorded, "`.",
      call. = FALSE
    )
  }
  species
}

.rc_medium_pattern <- function(name) {
  name <- tolower(trimws(as.character(name)))
  patterns <- c(
    glucose = "d[- ]?glucose|glucose|dextrose|(^|[^a-z0-9])glc([^a-z0-9]|$)",
    lactate = "l[- ]?lactate|lactate|lactic[ _-]?acid|(^|[^a-z0-9])lac([^a-z0-9]|$)",
    pyruvate = "pyruvate|pyruvic[ _-]?acid",
    acetate = "acetate|acetic[ _-]?acid",
    citrate = "citrate|citric[ _-]?acid",
    acetoacetate = "acetoacetate|acetoacetic[ _-]?acid",
    beta_hydroxybutyrate = "beta[- _]?hydroxybutyrate|3[- _]?hydroxybutyrate",
    alanine = "l[- ]?alanine|(^|[^a-z0-9])alanine([^a-z0-9]|$)",
    arginine = "l[- ]?arginine|(^|[^a-z0-9])arginine([^a-z0-9]|$)",
    asparagine = "l[- ]?asparagine|(^|[^a-z0-9])asparagine([^a-z0-9]|$)",
    aspartate = "l[- ]?aspartate|aspartic[ _-]?acid",
    cysteine = "l[- ]?cysteine|(^|[^a-z0-9])cysteine([^a-z0-9]|$)",
    cystine = "l[- ]?cystine|(^|[^a-z0-9])cystine([^a-z0-9]|$)",
    glutamate = "l[- ]?glutamate|glutamic[ _-]?acid",
    glutamine = "l[- ]?glutamine|(^|[^a-z0-9])glutamine([^a-z0-9]|$)",
    glycine = "(^|[^a-z0-9])glycine([^a-z0-9]|$)",
    histidine = "l[- ]?histidine|(^|[^a-z0-9])histidine([^a-z0-9]|$)",
    hydroxyproline = "hydroxy[- ]?l[- ]?proline|hydroxyproline",
    isoleucine = "l[- ]?isoleucine|(^|[^a-z0-9])isoleucine([^a-z0-9]|$)",
    leucine = "l[- ]?leucine|(^|[^a-z0-9])leucine([^a-z0-9]|$)",
    lysine = "l[- ]?lysine|(^|[^a-z0-9])lysine([^a-z0-9]|$)",
    methionine = "l[- ]?methionine|(^|[^a-z0-9])methionine([^a-z0-9]|$)",
    phenylalanine = "l[- ]?phenylalanine|(^|[^a-z0-9])phenylalanine([^a-z0-9]|$)",
    proline = "l[- ]?proline|(^|[^a-z0-9])proline([^a-z0-9]|$)",
    serine = "l[- ]?serine|(^|[^a-z0-9])serine([^a-z0-9]|$)",
    threonine = "l[- ]?threonine|(^|[^a-z0-9])threonine([^a-z0-9]|$)",
    tryptophan = "l[- ]?tryptophan|(^|[^a-z0-9])tryptophan([^a-z0-9]|$)",
    tyrosine = "l[- ]?tyrosine|(^|[^a-z0-9])tyrosine([^a-z0-9]|$)",
    valine = "l[- ]?valine|(^|[^a-z0-9])valine([^a-z0-9]|$)",
    biotin = "biotin|vitamin[ _-]?b7",
    choline = "choline",
    pantothenate = "pantothenate|pantothenic[ _-]?acid|vitamin[ _-]?b5",
    folate = "folate|folic[ _-]?acid|vitamin[ _-]?b9",
    niacinamide = "niacinamide|nicotinamide|vitamin[ _-]?b3",
    p_aminobenzoate = "p[- _]?aminobenzoate|para[- _]?aminobenzoic",
    pyridoxine = "pyridoxine|vitamin[ _-]?b6",
    riboflavin = "riboflavin|vitamin[ _-]?b2",
    thiamine = "thiamine|vitamin[ _-]?b1",
    vitamin_b12 = "cobalamin|vitamin[ _-]?b12",
    inositol = "myo[- _]?inositol|(^|[^a-z0-9])inositol([^a-z0-9]|$)",
    carnitine = "l[- ]?carnitine|(^|[^a-z0-9])carnitine([^a-z0-9]|$)",
    ethanolamine = "ethanolamine",
    creatine = "(^|[^a-z0-9])creatine([^a-z0-9]|$)",
    creatinine = "creatinine",
    urea = "(^|[^a-z0-9])urea([^a-z0-9]|$)",
    urate = "urate|uric[ _-]?acid",
    glutathione = "glutathione",
    oxygen = "oxygen|(^|[^a-z0-9])o2([^a-z0-9]|$)",
    carbon_dioxide = "carbon[ _-]?dioxide|(^|[^a-z0-9])co2([^a-z0-9]|$)",
    water = "(^|[^a-z0-9])water([^a-z0-9]|$)|(^|[^a-z0-9])h2o([^a-z0-9]|$)",
    bicarbonate = "bicarbonate|hydrogen[ _-]?carbonate",
    phosphate = "phosphate|orthophosphate|inorganic[ _-]?phosphate",
    sulfate = "sulfate|sulphate",
    chloride = "chloride",
    sodium = "(^|[^a-z0-9])sodium([^a-z0-9]|$)",
    potassium = "(^|[^a-z0-9])potassium([^a-z0-9]|$)",
    calcium = "(^|[^a-z0-9])calcium([^a-z0-9]|$)",
    magnesium = "(^|[^a-z0-9])magnesium([^a-z0-9]|$)",
    iron = "ferric|ferrous|(^|[^a-z0-9])iron([^a-z0-9]|$)"
  )
  if (name %in% names(patterns)) return(unname(patterns[[name]]))
  gsub("_", "[- _]?", name, fixed = TRUE)
}

.rc_medium_rows <- function(names, concentration_mM = NA_real_,
                            category = "nutrient", required = FALSE) {
  n <- length(names)
  concentration_mM <- rep(concentration_mM, length.out = n)
  category <- rep(category, length.out = n)
  required <- rep(required, length.out = n)
  data.frame(
    metabolite_name = names,
    metabolite_pattern = vapply(names, .rc_medium_pattern, character(1)),
    concentration_mM = as.numeric(concentration_mM),
    uptake_fraction = 1,
    category = category,
    target_exchange_flag = FALSE,
    required_match = as.logical(required),
    stringsAsFactors = FALSE
  )
}

.rc_medium_catalog <- function(preset_id, species = c("human", "mouse")) {
  species <- match.arg(species)
  amino_acids <- c(
    "alanine", "arginine", "asparagine", "aspartate", "cysteine",
    "cystine", "glutamate", "glutamine", "glycine", "histidine",
    "hydroxyproline", "isoleucine", "leucine", "lysine", "methionine",
    "phenylalanine", "proline", "serine", "threonine", "tryptophan",
    "tyrosine", "valine"
  )
  essential <- c(
    "histidine", "isoleucine", "leucine", "lysine", "methionine",
    "phenylalanine", "threonine", "tryptophan", "valine"
  )
  vitamins <- c(
    "biotin", "choline", "pantothenate", "folate", "niacinamide",
    "p_aminobenzoate", "pyridoxine", "riboflavin", "thiamine",
    "vitamin_b12", "inositol"
  )
  ions <- c(
    "sodium", "potassium", "calcium", "magnesium", "chloride",
    "bicarbonate", "phosphate", "sulfate", "iron"
  )
  physiologic_other <- c(
    "glucose", "lactate", "pyruvate", "acetate", "citrate",
    "acetoacetate", "beta_hydroxybutyrate", "carnitine", "ethanolamine",
    "creatine", "creatinine", "urea", "urate", "oxygen",
    "carbon_dioxide", "water"
  )
  physiologic <- rbind(
    .rc_medium_rows(
      amino_acids,
      category = "amino_acid",
      required = amino_acids %in% c(essential, "glutamine", "arginine")
    ),
    .rc_medium_rows(vitamins, category = "vitamin_or_cofactor"),
    .rc_medium_rows(ions, category = "inorganic_ion"),
    .rc_medium_rows(
      physiologic_other,
      category = "physiologic_metabolite",
      required = physiologic_other %in% c("glucose", "lactate", "oxygen")
    )
  )
  key_concentrations <- c(
    glucose = 5.0, lactate = 1.5, glutamine = 0.55, arginine = 0.11,
    alanine = 0.35, glycine = 0.25, serine = 0.12, pyruvate = 0.08,
    acetate = 0.10, citrate = 0.10, urate = 0.35, urea = 5.0,
    creatinine = 0.08, bicarbonate = 25, sodium = 140, potassium = 4.5,
    calcium = 1.2, magnesium = 0.8, chloride = 103, phosphate = 1.0
  )
  hit <- match(physiologic$metabolite_name, names(key_concentrations))
  physiologic$concentration_mM[!is.na(hit)] <- key_concentrations[hit[!is.na(hit)]]
  sensitivity_reference <- c(glucose = 25, lactate = 20, glutamine = 2)
  sensitivity_index <- match(
    physiologic$metabolite_name, names(sensitivity_reference)
  )
  sensitivity_rows <- !is.na(sensitivity_index)
  physiologic$uptake_fraction[sensitivity_rows] <- pmin(
    1,
    physiologic$concentration_mM[sensitivity_rows] /
      sensitivity_reference[sensitivity_index[sensitivity_rows]]
  )
  physiologic$target_exchange_flag[sensitivity_rows] <- TRUE

  if (preset_id %in% c("normal_human_plasma", "human_plasma", "hplm")) {
    return(physiologic)
  }
  if (identical(preset_id, "mouse_plasma")) {
    # Availability is based on the broad murine plasma/TIF metabolomics panel.
    # Concentrations remain provenance values rather than direct flux bounds.
    physiologic$concentration_mM[physiologic$metabolite_name == "glucose"] <- 7.5
    physiologic$concentration_mM[physiologic$metabolite_name == "lactate"] <- 2.0
    physiologic$uptake_fraction[sensitivity_rows] <- pmin(
      1,
      physiologic$concentration_mM[sensitivity_rows] /
        sensitivity_reference[sensitivity_index[sensitivity_rows]]
    )
    return(physiologic)
  }

  if (identical(preset_id, "rpmi1640")) {
    rpm_names <- c(
      "glycine", "hydroxyproline", "arginine", "asparagine", "aspartate",
      "cystine", "glutamate", "glutamine", "histidine", "isoleucine",
      "leucine", "lysine", "methionine", "phenylalanine", "proline",
      "serine", "threonine", "tryptophan", "tyrosine", "valine",
      "biotin", "choline", "pantothenate", "folate", "niacinamide",
      "p_aminobenzoate", "pyridoxine", "riboflavin", "thiamine",
      "vitamin_b12", "inositol", "calcium", "magnesium", "potassium",
      "sodium", "chloride", "bicarbonate", "phosphate", "glucose",
      "glutathione", "oxygen", "carbon_dioxide", "water"
    )
    rpm_conc <- c(
      0.13333334, 0.15267175, 1.1494253, 0.37878788, 0.15037593,
      0.20833333, 0.13605443, 2.0547945, 0.09677419, 0.3816794,
      0.3816794, 0.21857923, 0.10067114, 0.09090909, 0.17391305,
      0.2857143, 0.16806723, 0.024509804, 0.110497236, 0.17094018,
      8.1967213e-4, 0.021428572, 5.24109e-4, 0.002265, 0.008196721,
      0.00729927, 0.004854369, 5.319149e-4, 0.002967359,
      3.690037e-6, 0.19444445, 0.42372882, 0.40650406, 5.3333335,
      103.44827, 103.44827, 23.809525, 5.641791, 11.111111,
      0.003257329, NA, NA, NA
    )
    out <- .rc_medium_rows(
      rpm_names,
      rpm_conc,
      category = ifelse(
        rpm_names %in% amino_acids, "amino_acid",
        ifelse(rpm_names %in% vitamins, "vitamin_or_cofactor",
               ifelse(rpm_names %in% ions, "inorganic_ion", "other_component"))
      ),
      required = rpm_names %in% c("glucose", "glutamine", essential)
    )
    return(out)
  }

  if (identical(preset_id, "dmem_high_glucose")) {
    dmem_names <- c(
      "glycine", "arginine", "cystine", "glutamine", "histidine",
      "isoleucine", "leucine", "lysine", "methionine", "phenylalanine",
      "serine", "threonine", "tryptophan", "tyrosine", "valine",
      "choline", "pantothenate", "folate", "niacinamide", "pyridoxine",
      "riboflavin", "thiamine", "inositol", "calcium", "iron",
      "magnesium", "potassium", "sodium", "chloride", "bicarbonate",
      "phosphate", "glucose", "oxygen", "carbon_dioxide", "water"
    )
    dmem_conc <- c(
      0.4, 0.39810428, 0.20127796, 4.0, 0.2, 0.8015267, 0.8015267,
      0.7978142, 0.20134228, 0.4, 0.4, 0.79831934, 0.078431375,
      0.39846742, 0.8034188, 0.028571429, 0.008385744, 0.009070295,
      0.032786883, 0.019417476, 0.00106383, 0.011869436, 0.04,
      1.8018018, 2.4752476e-4, 0.8139166, 5.3333335, 110.344826,
      110.344826, 44.04762, 0.9057971, 25.0, NA, NA, NA
    )
    return(.rc_medium_rows(
      dmem_names,
      dmem_conc,
      category = ifelse(
        dmem_names %in% amino_acids, "amino_acid",
        ifelse(dmem_names %in% vitamins, "vitamin_or_cofactor",
               ifelse(dmem_names %in% ions, "inorganic_ion", "other_component"))
      ),
      required = dmem_names %in% c("glucose", "glutamine", essential)
    ))
  }

  base_id <- if (identical(species, "mouse")) "mouse_plasma" else
    "normal_human_plasma"
  if (preset_id %in% c(
    "high_glucose", "low_glucose", "high_lactate", "low_lactate",
    "low_glutamine"
  )) {
    out <- .rc_medium_catalog(base_id, species)
    target <- switch(
      preset_id,
      high_glucose = "glucose",
      low_glucose = "glucose",
      high_lactate = "lactate",
      low_lactate = "lactate",
      low_glutamine = "glutamine"
    )
    concentration <- switch(
      preset_id,
      high_glucose = 25,
      low_glucose = 1,
      high_lactate = 20,
      low_lactate = 0.5,
      low_glutamine = 0.05
    )
    reference_high <- switch(target, glucose = 25, lactate = 20, glutamine = 2)
    selected <- out$metabolite_name == target
    out$concentration_mM[selected] <- concentration
    out$uptake_fraction[selected] <- concentration / reference_high
    out$target_exchange_flag[selected] <- TRUE
    out$required_match[selected] <- TRUE
    return(out)
  }
  stop("Unsupported medium preset: ", preset_id, call. = FALSE)
}

.rc_medium_reference_catalog <- function() {
  data.frame(
    preset_id = c(
      "normal_human_plasma", "human_plasma", "hplm", "mouse_plasma",
      "rpmi1640", "dmem_high_glucose", "high_glucose", "low_glucose",
      "high_lactate", "low_lactate", "low_glutamine", "custom"
    ),
    species = c(
      rep("Homo sapiens", 3), "Mus musculus", rep("not species-specific", 2),
      rep("species-specific plasma background", 5), "user supplied"
    ),
    reference_label = c(
      rep("Cantor et al., Cell 2017 (HPLM); Psychogios et al., PLoS One 2011", 3),
      "Sullivan et al., eLife 2019; Wang et al., PNAS 2021",
      "Moore et al., JAMA 1967; Thermo Fisher RPMI-1640 formulation 11875",
      "Dulbecco and Freeman, Virology 1959; Thermo Fisher DMEM formulation 11965",
      "plasma background with explicit glucose sensitivity bound",
      "plasma background with explicit glucose sensitivity bound",
      "plasma background with explicit lactate sensitivity bound",
      "plasma background with explicit lactate sensitivity bound",
      "plasma background with explicit glutamine sensitivity bound",
      "user supplied"
    ),
    reference_doi = c(
      rep("10.1016/j.cell.2017.03.023;10.1371/journal.pone.0016957", 3),
      "10.7554/eLife.44235;10.1073/pnas.2102344118",
      "10.1001/jama.1967.03120080053007",
      "10.1016/0042-6822(59)90063-3",
      NA_character_, NA_character_, NA_character_, NA_character_, NA_character_,
      NA_character_
    ),
    evidence_scope = c(
      rep("adult human plasma-like polar nutrient availability", 3),
      "murine plasma and tumor interstitial-fluid nutrient availability",
      "serum-free RPMI-1640 basal formulation",
      "DMEM high-glucose basal formulation",
      "glucose sensitivity on a physiological plasma background",
      "glucose sensitivity on a physiological plasma background",
      "lactate sensitivity on a physiological plasma background",
      "lactate sensitivity on a physiological plasma background",
      "glutamine sensitivity on a physiological plasma background",
      "user-supplied extracellular environment"
    ),
    stringsAsFactors = FALSE
  )
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

.rc_compass_model_bound_medium <- function(gem, exchange_limit = 1) {
  validated <- rc_validate_gem(gem)
  if (!is.numeric(exchange_limit) || length(exchange_limit) != 1L ||
      !is.finite(exchange_limit) || exchange_limit <= 0) {
    stop("`exchange_limit` must be one positive finite number.", call. = FALSE)
  }
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem)
  }
  meta <- gem$reaction_meta[
    match(validated$reactions, as.character(gem$reaction_meta$reaction_id)),
    ,
    drop = FALSE
  ]
  exchange <- as.character(meta$reaction_id[as.character(meta$role) == "exchange"])
  exchange <- intersect(exchange, validated$reactions)
  if (!length(exchange)) {
    stop("No exchange reactions were identified in the GEM.", call. = FALSE)
  }
  index <- match(exchange, validated$reactions)
  original_lb <- as.numeric(validated$lb[index])
  original_ub <- as.numeric(validated$ub[index])
  data.frame(
    medium_scenario_id = "compass_model_bounds",
    exchange_reaction_id = exchange,
    metabolite_id = if ("metabolite_id" %in% colnames(meta)) {
      as.character(meta$metabolite_id[index])
    } else {
      NA_character_
    },
    preset_metabolite = NA_character_,
    concentration_mM = NA_real_,
    condition = "all",
    lb = pmax(original_lb, -exchange_limit),
    ub = pmin(original_ub, exchange_limit),
    available = TRUE,
    original_lb = original_lb,
    original_ub = original_ub,
    exchange_limit = exchange_limit,
    uptake_fraction = NA_real_,
    evidence_source = "gem_directionality_with_uniform_exchange_cap",
    assumption_level = "shared_model_defined_environment",
    target_exchange_flag = FALSE,
    concentration_used_for_rate_bound = FALSE,
    rate_bound_source = "original_gem_bounds_intersected_with_uniform_cap",
    stringsAsFactors = FALSE
  )
}

.rc_build_medium_preset <- function(
    gem, preset_id, species, exchange_limit, uptake_scale,
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
    ,
    drop = FALSE
  ]
  roles <- unique(trimws(as.character(exchange_roles)))
  roles <- roles[!is.na(roles) & nzchar(roles)]
  exchange_meta <- meta[as.character(meta$role) %in% roles, , drop = FALSE]
  if (!nrow(exchange_meta)) {
    stop("No exchange reactions found for the medium preset.", call. = FALSE)
  }
  if (is.null(compounds)) {
    compounds <- .rc_medium_catalog(preset_id, species)
  }
  required <- c(
    "metabolite_name", "metabolite_pattern", "concentration_mM",
    "uptake_fraction", "category", "target_exchange_flag", "required_match"
  )
  missing <- setdiff(required, colnames(compounds))
  if (length(missing)) {
    stop(
      "Preset compound table missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  annotation <- .rc_medium_annotation_text(exchange_meta)
  exchange_ids <- as.character(exchange_meta$reaction_id)
  matched_count <- integer(nrow(compounds))
  rows <- vector("list", nrow(compounds))
  scale <- .rc_medium_scale(uptake_scale, preset_id)
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
      as.numeric(compounds$uptake_fraction[[i]]) * scale
    )
    original_lb <- as.numeric(validated$lb[index])
    original_ub <- as.numeric(validated$ub[index])
    requested_lb <- -exchange_limit * fraction
    rows[[i]] <- data.frame(
      medium_scenario_id = preset_id,
      exchange_reaction_id = exchange_ids[hit],
      metabolite_id = if ("metabolite_id" %in% colnames(exchange_meta)) {
        as.character(exchange_meta$metabolite_id[hit])
      } else {
        NA_character_
      },
      preset_metabolite = as.character(compounds$metabolite_name[[i]]),
      nutrient_category = as.character(compounds$category[[i]]),
      concentration_mM = as.numeric(compounds$concentration_mM[[i]]),
      condition = as.character(condition %||% "all"),
      lb = pmax(original_lb, requested_lb),
      ub = pmin(original_ub, exchange_limit),
      available = TRUE,
      original_lb = original_lb,
      original_ub = original_ub,
      exchange_limit = exchange_limit,
      uptake_fraction = fraction,
      evidence_source = "literature_backed_medium_catalog",
      assumption_level = "availability_catalog_with_relative_uptake_cap",
      target_exchange_flag = as.logical(compounds$target_exchange_flag[[i]]),
      concentration_used_for_rate_bound =
        as.logical(compounds$target_exchange_flag[[i]]),
      rate_bound_source = if (isTRUE(compounds$target_exchange_flag[[i]])) {
        "relative_concentration_sensitivity_not_measured_flux"
      } else {
        "binary_availability_intersected_with_original_gem_directionality"
      },
      stringsAsFactors = FALSE
    )
  }
  unmatched <- compounds$metabolite_name[
    as.logical(compounds$required_match) & matched_count == 0L
  ]
  if (length(unmatched) && isTRUE(strict_preset_matching)) {
    stop(
      "Required medium components were not matched to GEM exchanges: ",
      paste(unmatched, collapse = ", "),
      ". Inspect GEM annotations or set `strict_preset_matching = FALSE` for a documented partial model.",
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
  reference_row <- reference[reference$preset_id == preset_id, , drop = FALSE]
  if (!nrow(reference_row) && preset_id %in% c(
    "high_glucose", "low_glucose", "high_lactate", "low_lactate",
    "low_glutamine"
  )) {
    reference_row <- reference[reference$preset_id == preset_id, , drop = FALSE]
  }
  if (!nrow(reference_row)) {
    reference_row <- data.frame(
      species = if (identical(species, "human")) "Homo sapiens" else "Mus musculus",
      reference_label = "user supplied",
      reference_doi = NA_character_,
      evidence_scope = "user-supplied extracellular environment",
      stringsAsFactors = FALSE
    )
  }
  output$species <- reference_row$species[[1L]]
  output$reference_label <- reference_row$reference_label[[1L]]
  output$reference_doi <- reference_row$reference_doi[[1L]]
  output$reference_pmid <- if ("reference_pmid" %in% colnames(reference_row)) {
    reference_row$reference_pmid[[1L]]
  } else {
    NA_character_
  }
  output$evidence_scope <- reference_row$evidence_scope[[1L]]
  rownames(output) <- NULL
  attr(output, "preset_diagnostics") <- data.frame(
    medium_scenario_id = preset_id,
    preset_metabolite = compounds$metabolite_name,
    nutrient_category = compounds$category,
    concentration_mM = compounds$concentration_mM,
    n_exchange_matches = matched_count,
    required_match = compounds$required_match,
    matched = matched_count > 0L,
    stringsAsFactors = FALSE
  )
  output
}

#' Build shared extracellular medium scenarios
#'
#' Named physiological and cell-culture media are complete availability
#' catalogs rather than small marker-metabolite lists. Concentrations are stored
#' as provenance and do not become fluxes. Bounds are always intersected with
#' the original GEM directionality, so a medium cannot open a blocked uptake or
#' secretion direction.
#'
#' @param gem A validated RegCompass GEM.
#' @param scenario One or more medium scenario identifiers. `"physiologic"`
#'   resolves to human or mouse plasma according to `species`.
#' @param species `"auto"`, `"human"`, or `"mouse"`.
#' @param custom_medium Exact reaction-level medium rows.
#' @param custom_metabolites Metabolite-pattern availability rows.
#' @param uptake_scale Non-negative global or named sensitivity multipliers.
#' @param condition_col Deprecated compatibility alias for `condition`.
#' @param exchange_roles Reaction roles treated as exchange reactions.
#' @param condition Shared condition label; canonical scoring requires `"all"`.
#' @param exchange_limit Absolute cap used for available exchange directions.
#' @param strict_preset_matching Stop when required medium components cannot be
#'   mapped to exchange annotations.
#' @return A medium constraint data frame.
#' @export
rc_make_medium_scenarios <- function(
    gem,
    scenario = "physiologic",
    species = c("auto", "human", "mouse"),
    custom_medium = NULL,
    custom_metabolites = NULL,
    uptake_scale = c(
      physiologic = 1,
      normal_human_plasma = 1,
      mouse_plasma = 1,
      rpmi1640 = 1,
      dmem_high_glucose = 1,
      high_glucose = 1,
      low_glucose = 1,
      high_lactate = 1,
      low_lactate = 1,
      low_glutamine = 1,
      permissive_all_exchange = 1,
      minimal = 0.1
    ),
    condition_col = NULL,
    exchange_roles = c("exchange"),
    condition = condition_col,
    exchange_limit = 1,
    strict_preset_matching = TRUE) {
  species <- .rc_infer_gem_species(gem, species)
  choices <- c(
    "physiologic", "compass_model_bounds", "permissive_all_exchange",
    "normal_human_plasma", "human_plasma", "hplm", "mouse_plasma",
    "rpmi1640", "dmem_high_glucose", "high_glucose", "low_glucose",
    "high_lactate", "low_lactate", "low_glutamine", "minimal",
    "blood_like", "culture_like", "tumor_low_glucose",
    "lactate_available", "custom"
  )
  scenario <- match.arg(scenario, choices = choices, several.ok = TRUE)
  aliases <- c(
    blood_like = if (identical(species, "human")) "normal_human_plasma" else
      "mouse_plasma",
    culture_like = "rpmi1640",
    tumor_low_glucose = "low_glucose",
    lactate_available = "high_lactate"
  )
  alias_hit <- intersect(scenario, names(aliases))
  if (length(alias_hit)) {
    warning(
      "Legacy medium aliases were mapped to canonical presets: ",
      paste(paste0(alias_hit, " -> ", aliases[alias_hit]), collapse = ", "),
      call. = FALSE
    )
    scenario <- unique(c(setdiff(scenario, alias_hit), unname(aliases[alias_hit])))
  }
  scenario[scenario == "physiologic"] <- if (identical(species, "human")) {
    "normal_human_plasma"
  } else {
    "mouse_plasma"
  }
  if (identical(species, "mouse") &&
      any(scenario %in% c("normal_human_plasma", "human_plasma", "hplm"))) {
    stop("Human plasma presets cannot be used with a Mouse-GEM run.", call. = FALSE)
  }
  if (identical(species, "human") && "mouse_plasma" %in% scenario) {
    stop("`mouse_plasma` requires a Mouse-GEM run.", call. = FALSE)
  }
  if (!is.null(custom_medium) && !is.null(custom_metabolites)) {
    stop("Supply only one of `custom_medium` or `custom_metabolites`.", call. = FALSE)
  }
  pieces <- list()
  if ("compass_model_bounds" %in% scenario) {
    pieces[[length(pieces) + 1L]] <- .rc_compass_model_bound_medium(
      gem,
      exchange_limit = exchange_limit
    )
  }
  if ("permissive_all_exchange" %in% scenario) {
    permissive <- .rc_compass_model_bound_medium(gem, exchange_limit)
    permissive$medium_scenario_id <- "permissive_all_exchange"
    permissive$evidence_source <- "technical_all_exchange_original_directions"
    permissive$assumption_level <- "technical_sensitivity_baseline"
    pieces[[length(pieces) + 1L]] <- permissive
  }
  preset_ids <- intersect(
    scenario,
    c(
      "normal_human_plasma", "human_plasma", "hplm", "mouse_plasma",
      "rpmi1640", "dmem_high_glucose", "high_glucose", "low_glucose",
      "high_lactate", "low_lactate", "low_glutamine"
    )
  )
  for (preset_id in preset_ids) {
    pieces[[length(pieces) + 1L]] <- .rc_build_medium_preset(
      gem = gem,
      preset_id = preset_id,
      species = species,
      exchange_limit = exchange_limit,
      uptake_scale = uptake_scale,
      condition = condition %||% "all",
      exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching
    )
  }
  if ("minimal" %in% scenario) {
    compounds <- .rc_medium_rows(
      c(
        "glucose", "glutamine", "arginine", "histidine", "isoleucine",
        "leucine", "lysine", "methionine", "phenylalanine", "threonine",
        "tryptophan", "valine", "oxygen", "water", "phosphate",
        "bicarbonate", "sodium", "potassium", "chloride"
      ),
      category = "minimal_required_nutrient",
      required = c(rep(TRUE, 13), rep(FALSE, 6))
    )
    pieces[[length(pieces) + 1L]] <- .rc_build_medium_preset(
      gem = gem,
      preset_id = "minimal",
      species = species,
      exchange_limit = exchange_limit,
      uptake_scale = uptake_scale,
      condition = condition %||% "all",
      exchange_roles = exchange_roles,
      strict_preset_matching = strict_preset_matching,
      compounds = compounds
    )
  }
  if ("custom" %in% scenario) {
    if (is.null(custom_medium) && is.null(custom_metabolites)) {
      stop(
        "`custom_medium` or `custom_metabolites` is required for `scenario = 'custom'`.",
        call. = FALSE
      )
    }
    if (!is.null(custom_medium)) {
      required <- c(
        "medium_scenario_id", "exchange_reaction_id", "lb", "ub", "available"
      )
      missing <- setdiff(required, colnames(custom_medium))
      if (length(missing)) {
        stop(
          "`custom_medium` missing columns: ", paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      custom <- custom_medium
      custom$exchange_reaction_id <- trimws(as.character(custom$exchange_reaction_id))
      custom$lb <- suppressWarnings(as.numeric(custom$lb))
      custom$ub <- suppressWarnings(as.numeric(custom$ub))
      custom$available <- as.logical(custom$available)
      if (anyNA(custom$exchange_reaction_id) ||
          any(!nzchar(custom$exchange_reaction_id)) ||
          any(!is.finite(custom$lb)) || any(!is.finite(custom$ub)) ||
          any(custom$lb > custom$ub) || anyNA(custom$available)) {
        stop(
          "Custom medium rows require valid reaction IDs, logical availability, and finite ordered bounds.",
          call. = FALSE
        )
      }
      optional_defaults <- list(
        metabolite_id = NA_character_, preset_metabolite = NA_character_,
        nutrient_category = "custom", concentration_mM = NA_real_,
        condition = as.character(condition %||% "all"),
        original_lb = NA_real_, original_ub = NA_real_,
        exchange_limit = exchange_limit, uptake_fraction = NA_real_,
        evidence_source = "user_supplied_custom_medium",
        assumption_level = "user_supplied",
        target_exchange_flag = FALSE,
        concentration_used_for_rate_bound = FALSE,
        rate_bound_source = "user_supplied_bound",
        species = if (identical(species, "human")) "Homo sapiens" else "Mus musculus",
        reference_label = "user supplied",
        reference_doi = NA_character_, reference_pmid = NA_character_,
        evidence_scope = "user-supplied reaction-level extracellular constraints"
      )
      for (name in names(optional_defaults)) {
        if (!name %in% colnames(custom)) custom[[name]] <- optional_defaults[[name]]
      }
      custom$condition[is.na(custom$condition) | !nzchar(custom$condition)] <-
        as.character(condition %||% "all")
      pieces[[length(pieces) + 1L]] <- custom
    } else {
      required <- c("metabolite_name", "metabolite_pattern", "available")
      missing <- setdiff(required, colnames(custom_metabolites))
      if (length(missing)) {
        stop(
          "`custom_metabolites` missing columns: ",
          paste(missing, collapse = ", "),
          call. = FALSE
        )
      }
      compounds <- custom_metabolites[custom_metabolites$available %in% TRUE, , drop = FALSE]
      if (!nrow(compounds)) {
        stop("`custom_metabolites` contains no available metabolites.", call. = FALSE)
      }
      defaults <- list(
        concentration_mM = NA_real_, uptake_fraction = 1,
        category = "custom", target_exchange_flag = FALSE,
        required_match = TRUE
      )
      for (name in names(defaults)) {
        if (!name %in% colnames(compounds)) compounds[[name]] <- defaults[[name]]
      }
      custom_reference <- data.frame(
        preset_id = "custom",
        species = if (identical(species, "human")) "Homo sapiens" else "Mus musculus",
        reference_label = as.character(compounds$reference_label[[1L]] %||% "user supplied"),
        reference_doi = as.character(compounds$reference_doi[[1L]] %||% NA_character_),
        reference_pmid = as.character(compounds$reference_pmid[[1L]] %||% NA_character_),
        evidence_scope = "user-supplied metabolite availability",
        stringsAsFactors = FALSE
      )
      pieces[[length(pieces) + 1L]] <- .rc_build_medium_preset(
        gem = gem,
        preset_id = "custom",
        species = species,
        exchange_limit = exchange_limit,
        uptake_scale = uptake_scale,
        condition = condition %||% "all",
        exchange_roles = exchange_roles,
        strict_preset_matching = strict_preset_matching,
        compounds = compounds,
        custom_reference = custom_reference
      )
    }
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
    "literature_catalog_with_original_gem_directionality_intersection"
  output
}

#' Apply medium constraints without expanding GEM directionality
#'
#' Every requested bound is intersected with the original GEM bounds. Closing
#' uptake therefore preserves an originally permitted secretion direction, and
#' no medium row can open a reaction direction that was blocked in the GEM.
#'
#' @param gem A validated GEM.
#' @param medium_table Medium rows from `rc_make_medium_scenarios()`.
#' @param condition Optional condition selector.
#' @param exchange_default_lb Default lower bound for unlisted exchanges.
#' @param exchange_default_ub Upper cap for unlisted exchanges.
#' @param allow_secretion Preserve originally permitted positive exchange flux.
#' @param strict Stop for unknown or non-exchange reaction IDs.
#' @return A list with the constrained GEM and bound diagnostics.
rc_apply_medium_constraints <- function(
    gem, medium_table, condition = NULL,
    exchange_default_lb = 0, exchange_default_ub = Inf,
    allow_secretion = TRUE, strict = TRUE) {
  if (!is.logical(allow_secretion) || length(allow_secretion) != 1L ||
      is.na(allow_secretion) ||
      !is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`allow_secretion` and `strict` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(exchange_default_lb) || length(exchange_default_lb) != 1L ||
      !is.finite(exchange_default_lb) ||
      !is.numeric(exchange_default_ub) || length(exchange_default_ub) != 1L ||
      is.na(exchange_default_ub) || exchange_default_lb > exchange_default_ub) {
    stop("Default exchange bounds must be ordered numeric scalars.", call. = FALSE)
  }
  validated <- rc_validate_gem(gem)
  reactions <- validated$reactions
  if (is.null(gem$reaction_meta) || !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem, medium_table = medium_table)
  }
  meta <- gem$reaction_meta[
    match(reactions, as.character(gem$reaction_meta$reaction_id)),
    ,
    drop = FALSE
  ]
  is_exchange <- as.character(meta$role) == "exchange"
  old_lb <- stats::setNames(as.numeric(validated$lb), reactions)
  old_ub <- stats::setNames(as.numeric(validated$ub), reactions)
  lb <- old_lb
  ub <- old_ub
  # Close unlisted uptake but preserve/cap only directions allowed by the GEM.
  lb[is_exchange] <- pmax(old_lb[is_exchange], exchange_default_lb)
  if (isTRUE(allow_secretion)) {
    ub[is_exchange] <- pmin(old_ub[is_exchange], exchange_default_ub)
  } else {
    ub[is_exchange] <- pmin(old_ub[is_exchange], 0, exchange_default_ub)
  }
  status <- stats::setNames(rep("not_exchange", length(reactions)), reactions)
  status[is_exchange] <- "exchange_default_uptake_closed"

  if (!is.null(medium_table)) {
    if (!is.data.frame(medium_table)) {
      stop("`medium_table` must be a data.frame.", call. = FALSE)
    }
    required <- c("exchange_reaction_id", "lb", "ub", "available")
    missing <- setdiff(required, colnames(medium_table))
    if (length(missing)) {
      stop(
        "`medium_table` missing columns: ", paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    medium <- medium_table
    medium$exchange_reaction_id <- trimws(as.character(medium$exchange_reaction_id))
    medium$condition <- if ("condition" %in% colnames(medium)) {
      as.character(medium$condition)
    } else {
      "all"
    }
    medium$condition[is.na(medium$condition) | !nzchar(medium$condition)] <- "all"
    keep <- medium$condition == "all"
    if (!is.null(condition)) keep <- keep | medium$condition == as.character(condition)
    medium <- medium[keep, , drop = FALSE]
    if (nrow(medium)) {
      medium$available <- as.logical(medium$available)
      medium$lb <- suppressWarnings(as.numeric(medium$lb))
      medium$ub <- suppressWarnings(as.numeric(medium$ub))
      if (anyNA(medium$exchange_reaction_id) ||
          any(!nzchar(medium$exchange_reaction_id)) ||
          anyNA(medium$available) ||
          any(!is.finite(medium$lb)) || any(!is.finite(medium$ub)) ||
          any(medium$lb > medium$ub)) {
        stop(
          "Medium rows require valid reaction IDs, logical availability, and finite ordered bounds.",
          call. = FALSE
        )
      }
      duplicate_key <- paste(
        medium$exchange_reaction_id,
        medium$condition,
        sep = "\001"
      )
      if (anyDuplicated(duplicate_key)) {
        stop(
          "`medium_table` contains duplicated reaction/condition rows.",
          call. = FALSE
        )
      }
      unknown <- setdiff(medium$exchange_reaction_id, reactions)
      if (length(unknown)) {
        message <- paste(
          "Medium exchange reactions missing from GEM:",
          paste(utils::head(unknown, 10L), collapse = ", ")
        )
        if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
      }
      medium <- medium[medium$exchange_reaction_id %in% reactions, , drop = FALSE]
      reaction_index <- match(medium$exchange_reaction_id, reactions)
      non_exchange <- medium$exchange_reaction_id[!is_exchange[reaction_index]]
      if (length(non_exchange)) {
        message <- paste(
          "Medium rows reference reactions not annotated as exchange:",
          paste(utils::head(unique(non_exchange), 10L), collapse = ", ")
        )
        if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
      }
      for (i in seq_len(nrow(medium))) {
        index <- reaction_index[[i]]
        if (!isTRUE(medium$available[[i]])) {
          lb[[index]] <- pmax(old_lb[[index]], exchange_default_lb)
          ub[[index]] <- if (isTRUE(allow_secretion)) {
            pmin(old_ub[[index]], exchange_default_ub)
          } else {
            pmin(old_ub[[index]], 0, exchange_default_ub)
          }
          status[[index]] <- "medium_unavailable_uptake_closed"
        } else {
          lb[[index]] <- max(old_lb[[index]], medium$lb[[i]])
          requested_ub <- if (isTRUE(allow_secretion)) {
            medium$ub[[i]]
          } else {
            min(medium$ub[[i]], 0)
          }
          ub[[index]] <- min(old_ub[[index]], requested_ub)
          status[[index]] <- "medium_available_intersection"
        }
      }
    }
  }
  if (any(lb > ub)) {
    bad <- reactions[lb > ub]
    stop(
      "Applied medium constraints produced lower bounds above upper bounds for: ",
      paste(utils::head(bad, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  gem$S <- validated$S
  gem$lb <- stats::setNames(as.numeric(lb), reactions)
  gem$ub <- stats::setNames(as.numeric(ub), reactions)
  gem$medium_policy <- "original_gem_directionality_intersection"
  diagnostics <- data.frame(
    reaction_id = reactions,
    old_lb = as.numeric(old_lb),
    old_ub = as.numeric(old_ub),
    new_lb = as.numeric(gem$lb),
    new_ub = as.numeric(gem$ub),
    lower_bound_expanded = gem$lb < old_lb,
    upper_bound_expanded = gem$ub > old_ub,
    medium_status = as.character(status),
    condition = condition %||% "all",
    stringsAsFactors = FALSE
  )
  if (any(diagnostics$lower_bound_expanded | diagnostics$upper_bound_expanded)) {
    stop("Medium application expanded the original GEM feasible region.", call. = FALSE)
  }
  list(gem = gem, medium_diagnostics = diagnostics)
}
