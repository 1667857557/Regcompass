.strict_fit_fixture <- function() {
  edge <- data.frame(
    edge_id = c("G||TF1||P1", "G||TF2||P2"),
    target = "G",
    tf = c("TF1", "TF2"),
    region = c("P1", "P2"),
    atac_feature_id = c("A1", "A2"),
    candidate_index = 1:2,
    stringsAsFactors = FALSE
  )
  pval <- c(0.001, 0.9, 0.02, 0.8)
  condition <- rep(c("A", "B"), each = 2L)
  padj <- unlist(lapply(split(pval, condition), stats::p.adjust,
                        method = "BH"), use.names = FALSE)
  estimate <- c(1, -0.2, -0.7, 0.1)
  significant <- padj < 0.05
  coefficient <- data.frame(
    edge_id = rep(edge$edge_id, 2L),
    target = "G",
    tf = rep(edge$tf, 2L),
    region = rep(edge$region, 2L),
    condition = condition,
    estimate = estimate,
    std_err = 0.1,
    statistic = estimate / 0.1,
    pval = pval,
    padj = padj,
    significant = significant,
    penalty_effect = ifelse(significant, estimate, 0),
    estimable = TRUE,
    zero_variance = FALSE,
    aliased = FALSE,
    direction = ifelse(estimate > 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )
  structure(list(
    schema_version = "pando_condition_grn_common_dictionary_v1",
    fit_engine = "two_stage_exact_edge_union_fixed_dictionary_glm",
    coefficient_scale = "shared_preprocessed_input_units_unscaled",
    inference_scope = "conditional_on_selected_edge_dictionary",
    cell_type = "T_cell",
    condition_levels = c("A", "B"),
    condition_col = "condition",
    cell_type_col = "cell_type",
    condition_cell_ids = list(A = c("a1", "a2"), B = c("b1", "b2")),
    edge_dictionary = edge,
    coefficients = coefficient,
    fit = data.frame(
      target = rep("G", 2), condition = c("A", "B"),
      rsq = c(0.8, 0.7), fit_status = "ok",
      stringsAsFactors = FALSE
    ),
    network_names = c(A = "net_A", B = "net_B"),
    padj_threshold = 0.05,
    adjust_method = "BH",
    scale = FALSE,
    interaction = ":",
    projection_effect_column = "penalty_effect",
    projection_policy = "padj_significant_effects_only"
  ), class = c("ConditionGRNFit", "list"))
}

test_that("strict validator requires the complete dictionary in every condition", {
  fit <- .strict_fit_fixture()
  expect_invisible(RegCompassR:::.rc_require_pando_condition_grn_fit(fit))
  incomplete <- fit
  incomplete$coefficients <- incomplete$coefficients[-1L, , drop = FALSE]
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(incomplete),
    "every frozen dictionary edge"
  )
})

test_that("significance and penalty effect are exactly estimable BH padj below 0.05", {
  fit <- .strict_fit_fixture()
  wrong_flag <- fit
  wrong_flag$coefficients$significant[[1L]] <- FALSE
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong_flag),
    "significant-edge flags"
  )
  wrong_effect <- fit
  wrong_effect$coefficients$penalty_effect[[1L]] <- 0
  expect_error(
    RegCompassR:::.rc_require_pando_condition_grn_fit(wrong_effect),
    "penalty_effect"
  )
})

test_that("condition contract is implemented once in functional source files", {
  root <- testthat::test_path("..", "..")
  expect_false(file.exists(file.path(
    root, "R", "zzzzzz_strict_condition_penalty.R"
  )))
  r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  source_text <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
  count_definition <- function(name) {
    pattern <- paste0("^", gsub("\\.", "\\\\.", name),
                      "[[:space:]]*<-[[:space:]]*function")
    sum(grepl(pattern, source_text))
  }
  expect_identical(count_definition(".rc_require_pando_condition_grn_fit"), 1L)
  expect_identical(count_definition(".rc_extract_condition_grn_contract"), 1L)
  expect_identical(count_definition(".rc_condition_pando_projection"), 1L)
  expect_false(any(grepl("strict_condition_penalty|_strict_base", source_text)))
})

test_that("package metadata pins the audited Pando head", {
  description <- read.dcf(testthat::test_path("..", "..", "DESCRIPTION"))
  remotes <- unname(description[1L, "Remotes"])
  expect_match(
    remotes,
    "Pando_regcompass@5c451e192831713ed955cd03276ddf47e89a4ec2",
    fixed = TRUE
  )
  expect_identical(
    unname(description[1L, "Config/RegCompass/PandoConditionEffectFilter"]),
    "BH-adjusted-p-below-0.05"
  )
})