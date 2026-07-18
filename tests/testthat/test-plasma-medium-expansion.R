make_exact_plasma_medium_gem <- function(species = "human") {
  metabolite_names <- c(
    "D-glucose", "hypoxanthine", "uridine", "taurine", "succinate",
    "glycerol", "O-acetylcarnitine", "2-oxoglutarate",
    "L-2-aminobutanoate", "N-acetylglycine",
    "ibuprofen-taurine conjugate"
  )
  metabolite_ids <- paste0("MAM_TEST_", seq_along(metabolite_names), "e")
  reaction_ids <- c(
    "EX_glucose", "EX_hypoxanthine", "EX_uridine", "EX_taurine",
    "EX_succinate", "EX_glycerol", "EX_acetylcarnitine",
    "EX_alpha_ketoglutarate", "EX_alpha_aminobutyrate",
    "EX_n_acetylglycine", "EX_taurine_conjugate"
  )
  S <- Matrix::Diagonal(length(reaction_ids), x = -1)
  dimnames(S) <- list(metabolite_ids, reaction_ids)
  list(
    S = S,
    lb = stats::setNames(rep(-1000, length(reaction_ids)), reaction_ids),
    ub = stats::setNames(rep(1000, length(reaction_ids)), reaction_ids),
    reaction_meta = data.frame(
      reaction_id = reaction_ids,
      reaction_name = c(
        paste("Exchange of", metabolite_names[-length(metabolite_names)]),
        "Exchange of taurine"
      ),
      role = "exchange",
      stringsAsFactors = FALSE
    ),
    metabolite_meta = data.frame(
      metabolite_id = metabolite_ids,
      name = metabolite_names,
      compartment = "e",
      stringsAsFactors = FALSE
    ),
    model_info = list(
      species = species,
      source = if (identical(species, "human")) {
        "SysBioChalmers/Human-GEM"
      } else {
        "SysBioChalmers/Mouse-GEM"
      }
    )
  )
}

test_that("plasma presets map literature nutrients to exact GEM metabolites", {
  medium <- rc_make_medium_scenarios(
    make_exact_plasma_medium_gem("human"),
    scenario = "normal_human_plasma",
    strict_preset_matching = FALSE
  )

  expected <- c(
    "glucose", "hypoxanthine", "uridine", "taurine", "succinate",
    "glycerol", "acetylcarnitine", "alpha_ketoglutarate",
    "alpha_aminobutyrate", "n_acetylglycine"
  )
  expect_true(all(expected %in% medium$preset_metabolite))
  expect_false("EX_taurine_conjugate" %in% medium$exchange_reaction_id)
  expect_true(all(
    medium$match_method[
      medium$preset_metabolite %in% expected
    ] == "exact_gem_metabolite_name"
  ))
  expect_true(all(
    grepl("e$", medium$metabolite_id[
      medium$preset_metabolite %in% expected
    ])
  ))
})

test_that("new HPLM concentrations are retained as provenance", {
  catalog <- .rc_medium_catalog("normal_human_plasma", "human")
  concentration <- stats::setNames(
    catalog$concentration_mM,
    catalog$metabolite_name
  )

  expect_equal(unname(concentration["hypoxanthine"]), 0.01)
  expect_equal(unname(concentration["uridine"]), 0.003001638)
  expect_equal(unname(concentration["taurine"]), 0.09000319)
  expect_equal(unname(concentration["succinate"]), 0.020001695)
  expect_equal(unname(concentration["glycerol"]), 0.11999132)
  expect_equal(unname(concentration["acetylcarnitine"]), 0.005002086)
})

test_that("sensitivity scenarios retain the expanded plasma background", {
  high <- rc_make_medium_scenarios(
    make_exact_plasma_medium_gem("human"),
    scenario = "high_glucose",
    strict_preset_matching = FALSE
  )

  expect_true(all(c(
    "hypoxanthine", "uridine", "taurine", "succinate"
  ) %in% high$preset_metabolite))
  expect_equal(
    high$concentration_mM[high$preset_metabolite == "glucose"],
    25
  )
})

test_that("mouse plasma uses detection evidence without borrowing human concentrations", {
  medium <- rc_make_medium_scenarios(
    make_exact_plasma_medium_gem("mouse"),
    scenario = "mouse_plasma",
    strict_preset_matching = FALSE
  )

  taurine <- medium[medium$preset_metabolite == "taurine", , drop = FALSE]
  expect_equal(nrow(taurine), 1)
  expect_true(is.na(taurine$concentration_mM))
  expect_equal(
    taurine$concentration_basis,
    "mouse_plasma_TIF_detection_without_shared_concentration"
  )
  expect_equal(taurine$component_reference_doi, "10.7554/eLife.44235")
})

test_that("preset diagnostics report exact aliases and unmatched nutrients", {
  medium <- rc_make_medium_scenarios(
    make_exact_plasma_medium_gem("human"),
    scenario = "normal_human_plasma",
    strict_preset_matching = FALSE
  )
  diagnostics <- attr(medium, "preset_diagnostics")

  expect_true(all(c(
    "gem_metabolite_aliases", "match_method", "concentration_basis",
    "component_reference_doi"
  ) %in% colnames(diagnostics)))
  expect_true(
    diagnostics$matched[diagnostics$preset_metabolite == "hypoxanthine"]
  )
  expect_false(
    diagnostics$matched[diagnostics$preset_metabolite == "folate"]
  )
})
