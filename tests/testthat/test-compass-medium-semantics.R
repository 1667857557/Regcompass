make_compass_medium_test_gem <- function() {
  reactions <- c(
    "EX_reverse", "EX_forward", "EX_unlisted",
    "EX_close_reverse", "EX_close_forward", "R_internal"
  )
  S <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 5, 6, 7),
    j = c(1, 2, 3, 4, 5, 6, 6),
    x = c(-1, 1, -1, -1, 1, -1, 1),
    dims = c(7, 6),
    dimnames = list(paste0("m", seq_len(7)), reactions)
  )
  list(
    S = S,
    lb = stats::setNames(rep(-1000, length(reactions)), reactions),
    ub = stats::setNames(rep(1000, length(reactions)), reactions),
    reaction_meta = data.frame(
      reaction_id = reactions,
      role = c(rep("exchange", 5), "internal"),
      stringsAsFactors = FALSE
    )
  )
}

test_that("COMPASS medium leaves omitted uptake capped instead of closed", {
  gem <- make_compass_medium_test_gem()
  medium <- data.frame(
    medium_scenario_id = "test",
    exchange_reaction_id = c(
      "EX_reverse", "EX_forward", "EX_close_reverse", "EX_close_forward"
    ),
    condition = "all",
    lb = c(-10, -1, 0, -1),
    ub = c(1, 4, 1, 0),
    available = c(TRUE, TRUE, FALSE, FALSE),
    exchange_limit = c(10, 4, 1, 1),
    uptake_fraction = 1,
    evidence_source = "literature_backed_medium_catalog",
    stringsAsFactors = FALSE
  )

  applied <- rc_apply_medium_constraints(gem, medium)
  expect_equal(unname(applied$gem$lb["EX_reverse"]), -10)
  expect_equal(unname(applied$gem$ub["EX_reverse"]), 1000)
  expect_equal(unname(applied$gem$lb["EX_forward"]), -1000)
  expect_equal(unname(applied$gem$ub["EX_forward"]), 4)
  expect_equal(unname(applied$gem$lb["EX_unlisted"]), -1)
  expect_equal(unname(applied$gem$ub["EX_unlisted"]), 1000)
  expect_equal(unname(applied$gem$lb["EX_close_reverse"]), 0)
  expect_equal(unname(applied$gem$ub["EX_close_reverse"]), 1000)
  expect_equal(unname(applied$gem$lb["EX_close_forward"]), -1000)
  expect_equal(unname(applied$gem$ub["EX_close_forward"]), 0)
  expect_equal(unname(applied$gem$lb["R_internal"]), -1000)
  expect_equal(unname(applied$gem$ub["R_internal"]), 1000)
  expect_identical(
    applied$gem$medium_semantics_version,
    "compass_exchange_bounds_v2"
  )

  contract <- .rc_assert_medium_bounds_only(gem, applied$gem, medium)
  expect_true(
    "EX_unlisted" %in% contract$changed_unlisted_exchange_reactions
  )
  expect_identical(contract$n_removed_reactions, 0L)
})

test_that("reaction-level custom medium bounds remain backward compatible", {
  gem <- make_compass_medium_test_gem()
  custom <- data.frame(
    medium_scenario_id = "custom",
    exchange_reaction_id = "EX_reverse",
    condition = "all",
    lb = -0.2,
    ub = 5,
    available = TRUE,
    stringsAsFactors = FALSE
  )
  applied <- rc_apply_medium_constraints(gem, custom)
  expect_equal(unname(applied$gem$lb["EX_reverse"]), -0.2)
  expect_equal(unname(applied$gem$ub["EX_reverse"]), 5)
})

test_that("COMPASS model-bound table caps uptake but preserves secretion", {
  gem <- make_compass_medium_test_gem()
  medium <- .rc_compass_model_bound_medium(gem, exchange_limit = 1)
  reverse <- medium[
    medium$exchange_reaction_id == "EX_reverse", , drop = FALSE
  ]
  forward <- medium[
    medium$exchange_reaction_id == "EX_forward", , drop = FALSE
  ]
  expect_equal(reverse$lb, -1)
  expect_equal(reverse$ub, 1000)
  expect_equal(forward$lb, -1000)
  expect_equal(forward$ub, 1)
  expect_true(all(medium$bound_scope == "uptake"))
})

test_that("medium audit allows exchange baseline changes only", {
  gem <- make_compass_medium_test_gem()
  medium <- data.frame(
    medium_scenario_id = "test",
    exchange_reaction_id = "EX_reverse",
    condition = "all",
    lb = -0.5,
    ub = 1000,
    available = TRUE,
    stringsAsFactors = FALSE
  )
  applied <- rc_apply_medium_constraints(gem, medium)
  expect_silent(.rc_assert_medium_bounds_only(gem, applied$gem, medium))
  invalid <- applied$gem
  invalid$ub[["R_internal"]] <- 5
  expect_error(
    .rc_assert_medium_bounds_only(gem, invalid, medium),
    "outside annotated exchange reactions"
  )
})
