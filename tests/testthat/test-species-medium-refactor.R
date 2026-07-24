make_species_medium_gem <- function(species = "human") {
  reactions <- c(
    "EX_glucose", "EX_lactate", "EX_glutamine", "EX_leucine",
    "EX_oxygen", "EX_blocked", "EX_uptake_only", "EX_secretion_only", "R1"
  )
  S <- Matrix::Matrix(
    matrix(c(rep(-1, 8), 1), nrow = 1,
           dimnames = list("m_e", reactions)),
    sparse = TRUE
  )
  list(
    S = S,
    lb = stats::setNames(
      c(-1000, -1000, -1000, -1000, -1000, 0, -1000, 0, 0),
      reactions
    ),
    ub = stats::setNames(
      c(1000, 1000, 1000, 1000, 1000, 0, 0, 5, 1000),
      reactions
    ),
    reaction_meta = data.frame(
      reaction_id = reactions,
      role = c(rep("exchange", 8), "internal"),
      metabolite_name = c(
        "D-glucose", "L-lactate", "L-glutamine", "L-leucine",
        "oxygen", "blocked nutrient", "uptake nutrient",
        "secretion nutrient", NA
      ),
      stringsAsFactors = FALSE
    ),
    model_info = list(
      species = species,
      source = if (identical(species, "human")) {
        "SysBioChalmers/Human-GEM"
      } else {
        "SysBioChalmers/Mouse-GEM"
      },
      version = if (identical(species, "human")) "2.0.0" else "1.8.0"
    )
  )
}

test_that("species model routing is explicit and pinned", {
  human <- .rc_species_gem_spec("human")
  mouse <- .rc_species_gem_spec("mouse")
  expect_equal(human$source, "SysBioChalmers/Human-GEM")
  expect_equal(human$version, "2.0.0")
  expect_equal(human$taxonomy_id, "9606")
  expect_equal(mouse$source, "SysBioChalmers/Mouse-GEM")
  expect_equal(mouse$version, "1.8.0")
  expect_equal(mouse$taxonomy_id, "10090")
  expect_identical(eval(formals(rc_prepare_gem)$species), c("human", "mouse"))
  expect_identical(eval(formals(rc_prepare_human2_gem)$version), "2.0.0")
  expect_identical(eval(formals(rc_prepare_mouse_gem)$version), "1.8.0")
})

test_that("Mouse-GEM GPR rules retain mouse symbols directly", {
  root <- tempfile("mouse-gem-fixture-")
  dir.create(file.path(root, "model"), recursive = TRUE)
  writeLines(
    c(
      "- reactions:",
      "    - id: 'MAR00001'",
      "      gene_reaction_rule: 'Adh1 or (Adh4 and Aldh2)'"
    ),
    file.path(root, "model", "Mouse-GEM.yml")
  )
  utils::write.table(
    data.frame(rxns = "MAR00001"),
    file.path(root, "model", "reactions.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  prepared <- rc_prepare_species_gpr_table(
    root, species = "mouse", gene_format = "symbol"
  )
  groups <- split(prepared$gpr_table$gene, prepared$gpr_table$and_group_id)
  keys <- sort(vapply(groups, function(x) paste(sort(x), collapse = "+"),
                      character(1)))
  expect_equal(unname(keys), c("Adh1", "Adh4+Aldh2"))
  expect_false(any(grepl("ENSG", prepared$gpr_table$gene)))
})

test_that("literature media are complete catalogs rather than marker lists", {
  human_plasma <- .rc_medium_catalog("normal_human_plasma", "human")
  mouse_plasma <- .rc_medium_catalog("mouse_plasma", "mouse")
  rpm <- .rc_medium_catalog("rpmi1640", "human")
  dmem <- .rc_medium_catalog("dmem_high_glucose", "human")
  expect_gte(nrow(human_plasma), 55)
  expect_gte(nrow(mouse_plasma), 55)
  expect_gte(nrow(rpm), 40)
  expect_gte(nrow(dmem), 30)
  expect_true(all(c(
    "glucose", "glutamine", "leucine", "methionine", "tryptophan",
    "choline", "folate", "thiamine", "sodium", "potassium",
    "chloride", "bicarbonate", "phosphate", "oxygen", "water"
  ) %in% rpm$metabolite_name))
  expect_true(all(c(
    "glucose", "glutamine", "leucine", "methionine", "tryptophan",
    "choline", "folate", "thiamine", "sodium", "calcium",
    "magnesium", "chloride", "bicarbonate", "phosphate", "oxygen",
    "water"
  ) %in% dmem$metabolite_name))
})

test_that("physiological default follows GEM species", {
  human <- rc_make_medium_scenarios(
    make_species_medium_gem("human"), strict_preset_matching = FALSE
  )
  mouse <- rc_make_medium_scenarios(
    make_species_medium_gem("mouse"), strict_preset_matching = FALSE
  )
  expect_true(all(human$medium_scenario_id == "normal_human_plasma"))
  expect_true(all(mouse$medium_scenario_id == "mouse_plasma"))
  expect_identical(attr(human, "species"), "human")
  expect_identical(attr(mouse, "species"), "mouse")
})

test_that("medium application never expands blocked GEM directions", {
  gem <- make_species_medium_gem("human")
  medium <- data.frame(
    medium_scenario_id = "test",
    exchange_reaction_id = c(
      "EX_blocked", "EX_uptake_only", "EX_secretion_only", "EX_glucose"
    ),
    condition = "all",
    lb = c(-1, -1, -1, -0.2),
    ub = c(1, 1, 10, 1),
    available = TRUE,
    stringsAsFactors = FALSE
  )
  constrained <- rc_apply_medium_constraints(gem, medium)
  expect_equal(unname(constrained$gem$lb["EX_blocked"]), 0)
  expect_equal(unname(constrained$gem$ub["EX_blocked"]), 0)
  expect_equal(unname(constrained$gem$ub["EX_uptake_only"]), 0)
  expect_equal(unname(constrained$gem$lb["EX_secretion_only"]), 0)
  expect_equal(unname(constrained$gem$ub["EX_secretion_only"]), 5)
  expect_equal(unname(constrained$gem$lb["EX_glucose"]), -0.2)
  expect_false(any(
    constrained$medium_diagnostics$lower_bound_expanded |
      constrained$medium_diagnostics$upper_bound_expanded
  ))
})

test_that("canonical Layer 2 owns a persistent model cache", {
  workflow_text <- paste(deparse(body(rc_regcompass_step_layer2)), collapse = "\n")
  micro_text <- paste(deparse(body(rc_run_microcompass)), collapse = "\n")
  expect_match(workflow_text, 'file.path(\n    outdir, "model_cache", model_mode', fixed = TRUE)
  expect_match(workflow_text, 'layer2_args$model_params$cache_dir', fixed = TRUE)
  expect_match(micro_text, "model_file_manifest.rds", fixed = TRUE)
  expect_match(micro_text, "tools::md5sum", fixed = TRUE)
})

test_that("one-shot species argument routes setup by species", {
  expect_identical(
    eval(formals(rc_run_regcompass_one_shot)$species), c("human", "mouse")
  )
  body_text <- paste(deparse(body(rc_run_regcompass_one_shot)), collapse = "\n")
  expect_match(body_text, "rc_prepare_gem", fixed = TRUE)
  expect_match(body_text, "species = species", fixed = TRUE)
  expect_match(body_text, "rc_make_medium_scenarios", fixed = TRUE)
  expect_match(body_text, "Mouse-GEM|1.8.0")
})
