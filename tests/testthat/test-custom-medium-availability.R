test_that("unavailable custom metabolites remain unavailable after GEM mapping", {
  reactions <- c("EX_reverse", "EX_forward", "R1")
  S <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(1, 2, 3, 3),
    x = c(-1, 1, -1, 1),
    dims = c(4, 3),
    dimnames = list(paste0("m", seq_len(4)), reactions)
  )
  gem <- list(
    S = S,
    lb = stats::setNames(c(-1000, -1000, 0), reactions),
    ub = stats::setNames(c(1000, 1000, 1000), reactions),
    model_info = list(species = "human"),
    reaction_meta = data.frame(
      reaction_id = reactions,
      role = c("exchange", "exchange", "internal"),
      metabolite_name = c("reversefuel", "forwardfuel", NA_character_),
      stringsAsFactors = FALSE
    )
  )

  custom <- data.frame(
    metabolite_name = c("reversefuel", "forwardfuel"),
    available = c(FALSE, FALSE),
    uptake_fraction = c(1, 1),
    target_exchange_flag = TRUE,
    required_match = TRUE,
    stringsAsFactors = FALSE
  )
  medium <- rc_make_medium_scenarios(
    gem,
    scenario = NULL,
    species = "human",
    custom_metabolites = custom,
    exchange_limit = 2
  )

  expect_true(all(medium$available %in% FALSE))
  expect_true(all(medium$uptake_fraction == 0))

  applied <- rc_apply_medium_constraints(gem, medium)
  expect_equal(unname(applied$gem$lb["EX_reverse"]), 0)
  expect_equal(unname(applied$gem$ub["EX_forward"]), 0)
})
