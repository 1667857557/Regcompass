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

test_that("human plasma maps authoritative HPLM nutrients exactly", {
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
  expect_true(all(
    medium$medium_background_id == "authoritative_HPLM_2017_2021"
  ))
})

test_that("updated HPLM components cite Cell Metabolism 2021", {
  medium <- rc_make_medium_scenarios(
    make_exact_plasma_medium_gem("human"),
    scenario = "normal_human_plasma",
    strict_preset_matching = FALSE
  )
  updated <- medium$preset_metabolite %in%
    c("uridine", "acetylcarnitine", "alpha_ketoglutarate")
  expect_true(all(
    medium$component_reference_doi[updated] == "10.1016/j.cmet.2021.02.005"
  ))
  expect_true(all(grepl(
    "updated_HPLM", medium$concentration_basis[updated], fixed = TRUE
  )))
})

test_that("challenge backgrounds retain authoritative HPLM small molecules", {
  high <- rc_make_medium_scenarios(
    make_exact_plasma_medium_gem("human"),
    scenario = "high_glucose",
    strict_preset_matching = FALSE
  )

  expect_equal(
    high$concentration_mM[high$preset_metabolite == "glucose"],
    25
  )
  expect_true(all(
    high$medium_background_id == "authoritative_HPLM_2017_2021"
  ))
  expect_true(all(
    high$scenario_construction ==
      "authoritative_HPLM_background_plus_named_nutrient_override"
  ))
  expect_true(all(c(
    "hypoxanthine", "uridine", "taurine", "succinate"
  ) %in% high$preset_metabolite))
})

test_that("mouse plasma exposes only the conservative Nature-supported subset", {
  medium <- rc_make_medium_scenarios(
    make_exact_plasma_medium_gem("mouse"),
    scenario = "mouse_plasma",
    strict_preset_matching = FALSE
  )

  expect_true(all(c(
    "glucose", "hypoxanthine", "uridine"
  ) %in% medium$preset_metabolite))
  expect_false(any(c(
    "taurine", "succinate", "glycerol", "acetylcarnitine",
    "alpha_ketoglutarate", "alpha_aminobutyrate", "n_acetylglycine"
  ) %in% medium$preset_metabolite))
  expect_true(all(
    medium$medium_background_id == "Abbott_2026_Nature_mouse_plasma"
  ))
  expect_true(all(grepl(
    "10.1038/s41586-025-09898-9",
    medium$reference_doi,
    fixed = TRUE
  )))
  availability <- medium[medium$preset_metabolite %in%
                           c("hypoxanthine", "uridine"), , drop = FALSE]
  expect_true(all(is.na(availability$concentration_mM)))
  expect_true(all(
    availability$component_reference_doi == "10.1038/s41586-025-09898-9"
  ))
  expect_false(any(availability$target_exchange_flag))
})

test_that("preset diagnostics retain source-specific component provenance", {
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
  expect_true(
    diagnostics$matched[diagnostics$preset_metabolite == "uridine"]
  )
})
