test_that("shared column helper has one stable calling convention", {
  x <- data.frame(a = 1, b = 2)
  expect_equal(.rc_first_existing_col(x, c("z", "b")), "b")
  expect_equal(.rc_first_existing_col(x, "z", fallback = "a"), "a")
  expect_null(.rc_first_existing_col(x, "z"))
})

test_that("microCOMPASS row-id parser has one implementation", {
  r_dir <- testthat::test_path("../../R")
  sources <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  definitions <- unlist(lapply(sources, function(path) {
    grep(
      "^rc_parse_microcompass_row_id[[:space:]]*<-[[:space:]]*function",
      readLines(path, warn = FALSE),
      value = TRUE
    )
  }), use.names = FALSE)
  expect_length(definitions, 1L)
})

test_that("removed compatibility arguments stay removed", {
  expect_false("min_direct" %in% names(formals(rc_q95_calibrate)))
  expect_false("min_direct" %in% names(formals(rc_run_layer1_capacity)))
  expect_false("require_model_info" %in% names(formals(rc_read_gem)))
  expect_false("require_model_info" %in% names(formals(rc_prepare_human2_gem)))
  expect_false("require_model_info" %in% names(formals(rc_prepare_human2_gem_v12)))
})

test_that("hard-min GPR capacity requires every AND subunit", {
  gpr <- list(R1 = list(c("g1", "g2")))
  score <- matrix(1, nrow = 1, dimnames = list("g1", "u1"))
  out <- rc_hard_min_capacity(gpr, score)
  expect_true(is.na(out["R1", "u1"]))
})

test_that("GPR parsing rejects duplicate reactions and trims genes", {
  expect_error(
    rc_parse_gpr_table(data.frame(
      reaction_id = c("R1", "R1"), gpr = c("g1", "g2")
    )),
    "one row per reaction"
  )
  parsed <- rc_parse_gpr_table(data.frame(
    reaction_id = " R1 ", and_group_id = " 1 ", gene = " G1 "
  ))
  expect_identical(names(parsed), "R1")
  expect_identical(parsed[[1]][[1]], "g1")
})

test_that("GEM validation enforces zero-column and bound contracts", {
  S <- matrix(0, nrow = 1, ncol = 1, dimnames = list("m", "R"))
  expect_error(
    rc_validate_gem(list(S = S, lb = c(R = 0), ub = c(R = 1)),
                    allow_zero_support = FALSE),
    "Zero-stoichiometry"
  )
  expect_error(
    rc_make_gem(S, lb = c(UNKNOWN = 0), ub = c(R = 1)),
    "Unknown reaction IDs"
  )
})

test_that("reaction role inference switches act independently", {
  S <- matrix(1, nrow = 1, dimnames = list("m_c", "EX_fake"))
  gem <- rc_make_gem(S, lb = c(EX_fake = 0), ub = c(EX_fake = 1))
  no_id <- rc_annotate_reaction_roles(
    gem, infer_from_id = FALSE, infer_from_stoichiometry = FALSE,
    infer_from_compartment = FALSE, overwrite_existing = TRUE
  )
  expect_equal(no_id$reaction_meta$role, "unknown")
  with_id <- rc_annotate_reaction_roles(
    gem, infer_from_id = TRUE, infer_from_stoichiometry = FALSE,
    infer_from_compartment = FALSE, overwrite_existing = TRUE
  )
  expect_equal(with_id$reaction_meta$role, "exchange")
})

test_that("medium constraints cannot reopen secretion when disabled", {
  S <- matrix(1, nrow = 1, dimnames = list("m_e", "EX_m"))
  gem <- rc_make_gem(
    S, lb = c(EX_m = -1000), ub = c(EX_m = 1000),
    reaction_meta = data.frame(
      reaction_id = "EX_m", role = "exchange", role_source = "curated"
    )
  )
  medium <- data.frame(
    exchange_reaction_id = "EX_m", lb = -10, ub = 1000,
    available = TRUE, condition = "all"
  )
  out <- rc_apply_medium_constraints(gem, medium, allow_secretion = FALSE)
  expect_equal(out$gem$ub[["EX_m"]], 0)
})

test_that("Q95 all-missing flag is global rather than per stratum", {
  C <- matrix(
    c(1, NA, NA, NA), nrow = 1,
    dimnames = list("R1", c("u1", "u2", "u3", "u4"))
  )
  meta <- data.frame(
    pool_id = colnames(C), cell_type = c("A", "A", "B", "B")
  )
  out <- suppressWarnings(rc_q95_shrink(C, meta, "cell_type"))
  expect_false(any(out$Q$all_missing_reaction_flag))
  expect_true(out$Q$stratum_missing_reaction_flag[out$Q$stratum == "B"])
})

test_that("partial-group confidence and threshold parameters are effective", {
  gpr <- list(R1 = list(c("g1", "g2")))
  confidence <- matrix(0.2, nrow = 1, dimnames = list("g1", "u1"))
  complete <- rc_reaction_confidence_gpr_aware(
    gpr, confidence, missing_group_policy = "complete_group"
  )
  partial <- rc_reaction_confidence_gpr_aware(
    gpr, confidence, missing_group_policy = "partial_group",
    low_confidence_threshold = 0.5
  )
  expect_true(is.na(complete$reaction_confidence))
  expect_equal(partial$reaction_confidence, 0.2)
  expect_true(partial$low_confidence_reaction_flag)
  expect_equal(partial$n_and_groups_eligible, 1L)
})

test_that("metacell metadata rejects mixed biological strata", {
  membership <- data.frame(
    metacell_id = c("M1", "M1"), cell_id = c("c1", "c2"),
    sample_id = c("S1", "S2"), condition = "ctrl", cell_type = "T"
  )
  expect_error(rc_build_metacell_metadata(membership), "mixes metadata")
})

test_that("pseudobulk rejects duplicated cell assignments", {
  counts <- Matrix::Matrix(
    matrix(c(1, 2), nrow = 1, dimnames = list("g", c("c1", "c2"))),
    sparse = TRUE
  )
  map <- data.frame(
    pool_id = c("p1", "p2"), cell_id = c("c1", "c1")
  )
  expect_error(rc_unit_bulk_counts(counts, map), "only one active pool")
})

test_that("GRN hard core requires one complete GPR AND group", {
  nodes <- data.frame(
    sample_id = "S1", module_id = "M1", gene = c("G1", "G3")
  )
  gpr <- data.frame(
    reaction_id = c("R1", "R1", "R2"),
    and_group_id = c(1, 1, 1), gene = c("G1", "G2", "G3")
  )
  mapped <- rc_map_meta_module_core_reactions(nodes, gpr)
  expect_false(unique(mapped$is_core[mapped$reaction_id == "R1"]))
  expect_true(unique(mapped$is_core[mapped$reaction_id == "R2"]))
})

test_that("meta-module expansion ignores incomplete GPR mappings", {
  S <- matrix(
    c(1, 0, 0, 1), nrow = 2,
    dimnames = list(c("m1", "m2"), c("R1", "R2"))
  )
  gem <- rc_make_gem(
    S,
    lb = c(R1 = 0, R2 = 0),
    ub = c(R1 = 1000, R2 = 1000),
    reaction_meta = data.frame(
      reaction_id = c("R1", "R2"),
      subsystem = c("A", "B"),
      stringsAsFactors = FALSE
    )
  )
  core <- data.frame(
    sample_id = "S1",
    module_id = "M1",
    gene = c("G1", "G3"),
    reaction_id = c("R1", "R2"),
    is_core = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  expanded <- rc_expand_meta_module_reactions(gem, core)
  expect_false("R1" %in% expanded$reaction_membership$reaction_id)
  expect_true("R2" %in% expanded$reaction_membership$reaction_id)
})

test_that("Layer 2 reaction helpers respect role and mixed GPR evidence", {
  meta <- data.frame(
    reaction_id = c("R1", "R2"), role = c("transport", "internal"),
    gpr = c("g1", "")
  )
  expect_equal(rc_layer2_reaction_type(meta), c("transport", "other"))
  C <- matrix(NA_real_, nrow = 2, ncol = 1,
              dimnames = list(c("R1", "R2"), "u1"))
  F <- C
  C["R2", "u1"] <- 1
  expect_equal(rc_layer2_has_gpr(meta, C, F), c(TRUE, TRUE))
})

test_that("microCOMPASS differential testing uses biological samples", {
  unit_meta <- data.frame(
    unit_id = paste0("u", 1:7),
    sample_id = c("A1", "A1", "A1", "A2", "B1", "B2", "B2"),
    condition = c("A", "A", "A", "A", "B", "B", "B"),
    cell_type = "T"
  )
  # Sample medians: A1=0, A2=2, B1=4, B2=6; B-A effect is 4.
  score <- matrix(
    c(0, 0, 100, 2, 4, 6, 6), nrow = 1,
    dimnames = list("reaction=R1::direction=forward::medium=base", unit_meta$unit_id)
  )
  out <- rc_test_microcompass_differential(
    list(score = score, unit_meta = unit_meta),
    method = "lm", min_samples_per_group = 2,
    preferred_min_samples_per_group = 2,
    strict_replicate_design = TRUE
  )
  expect_equal(out$effect_size, 4, tolerance = 1e-8)
  expect_equal(out$n_samples_per_group, "A=2;B=2")
})

test_that("Layer 2 linear models fit one row per biological sample", {
  score <- matrix(
    c(0, 10, 2, 4, 6), nrow = 1,
    dimnames = list("R1", paste0("u", 1:5))
  )
  meta <- data.frame(
    unit_id = colnames(score),
    sample_id = c("S1", "S1", "S2", "S3", "S4"),
    condition = c("A", "A", "A", "B", "B"),
    cell_type = "T"
  )
  fits <- rc_layer2_lm_by_reaction(
    score, meta, L2_score ~ condition
  )
  expect_equal(stats::nobs(fits[["T::R1"]]), 4L)
})

test_that("blocked hard core is explicitly marked as no allowed direction", {
  skip_if_not(requireNamespace("highs", quietly = TRUE))
  S <- matrix(1, nrow = 1, dimnames = list("m", "R"))
  gem <- rc_make_gem(
    S, lb = c(R = 0), ub = c(R = 0),
    reaction_meta = data.frame(
      reaction_id = "R", role = "internal", role_source = "curated"
    )
  )
  membership <- data.frame(
    sample_id = "S1", module_id = "M1", reaction_id = "R", is_core = TRUE
  )
  model <- rc_build_meta_module_gem(
    gem, membership, membership, "S1", "M1", solver = "highs"
  )
  expect_equal(model$target_status, "no_allowed_direction")
  expect_equal(model$closure_diagnostics$completion_status, "no_allowed_direction")
})

test_that("condition-specific full-GEM caches are distinct", {
  S <- matrix(1, nrow = 1, dimnames = list("m_e", "EX_m"))
  gem <- rc_make_gem(
    S, lb = c(EX_m = -1000), ub = c(EX_m = 1000),
    reaction_meta = data.frame(
      reaction_id = "EX_m", role = "exchange", role_source = "curated"
    )
  )
  medium <- data.frame(
    medium_scenario_id = "custom", exchange_reaction_id = "EX_m",
    lb = c(-1, -10), ub = 1000, available = TRUE,
    condition = c("A", "B")
  )
  dirs <- data.frame(reaction_id = "EX_m", target_direction = "reverse")
  cache <- rc_build_full_gem_cache(
    gem, dirs, medium, conditions = c("A", "B")
  )
  summary <- attr(cache, "summary")
  expect_setequal(summary$condition, c("A", "B"))
  models <- lapply(summary$file, readRDS)
  expect_setequal(vapply(models, function(x) x$lb[["EX_m"]], numeric(1)), c(-1, -10))
})
