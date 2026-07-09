test_that("Human-GEM metabolite compartment is parsed from metabolite IDs", {
  ids <- c("MAM00001c", "MAM00002e", "MAM00003m")
  expect_equal(.rc_humangem_met_compartment_from_id(ids), c("c", "e", "m"))
})

test_that("single extracellular metabolite reactions are exchange in Human-GEM style IDs", {
  S <- matrix(c(1, 0, 0, 1), nrow = 2,
              dimnames = list(c("MAM00001e", "MAM00002c"), c("MAR_EXLIKE", "MAR_INTERNAL")))
  met <- data.frame(metabolite_id = rownames(S), compartment = c("e", "c"), stringsAsFactors = FALSE)
  gem <- rc_make_gem(S, lb = c(MAR_EXLIKE = -1000, MAR_INTERNAL = 0),
                     ub = c(MAR_EXLIKE = 1000, MAR_INTERNAL = 1000), metabolite_meta = met)
  gem <- rc_annotate_reaction_roles(gem, overwrite_existing = TRUE)
  roles <- stats::setNames(gem$reaction_roles$role, gem$reaction_roles$reaction_id)
  expect_equal(roles[["MAR_EXLIKE"]], "exchange")
})

test_that("cross-compartment reactions are transport", {
  S <- matrix(c(-1, 1), nrow = 2, dimnames = list(c("MAM00001c", "MAM00001e"), "MAR_TRANSPORT"))
  met <- data.frame(metabolite_id = rownames(S), compartment = c("c", "e"), stringsAsFactors = FALSE)
  gem <- rc_make_gem(S, lb = c(MAR_TRANSPORT = 0), ub = c(MAR_TRANSPORT = 1000), metabolite_meta = met)
  gem <- rc_annotate_reaction_roles(gem, overwrite_existing = TRUE)
  roles <- stats::setNames(gem$reaction_roles$role, gem$reaction_roles$reaction_id)
  expect_equal(roles[["MAR_TRANSPORT"]], "transport")
})
