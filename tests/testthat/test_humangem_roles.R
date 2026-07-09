test_that("Human-GEM metabolite compartment is parsed from metabolite IDs", {
  ids <- c("MAM00001c", "MAM00002e", "MAM00003m", "MAM00004[e]")
  expect_equal(.rc_humangem_met_compartment_from_id(ids), c("c", "e", "m", "e"))
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

test_that("extracellular transport equations are not mislabeled as exchange", {
  S <- matrix(c(-1, 1), nrow = 2, dimnames = list(c("MAM00001c", "MAM00001e"), "MAR_EXT_TRANSPORT"))
  met <- data.frame(metabolite_id = rownames(S), compartment = c("c", "e"), stringsAsFactors = FALSE)
  meta <- data.frame(reaction_id = "MAR_EXT_TRANSPORT", equation = "MAM00001[c] <=> MAM00001[e]", stringsAsFactors = FALSE)
  gem <- rc_make_gem(S, lb = c(MAR_EXT_TRANSPORT = -1000), ub = c(MAR_EXT_TRANSPORT = 1000), reaction_meta = meta, metabolite_meta = met)
  gem <- rc_annotate_reaction_roles(gem, overwrite_existing = TRUE)
  roles <- stats::setNames(gem$reaction_roles$role, gem$reaction_roles$reaction_id)
  expect_equal(roles[["MAR_EXT_TRANSPORT"]], "transport")
})

test_that("medium table can identify exchange reactions with Human-GEM style IDs", {
  S <- matrix(1, nrow = 1, dimnames = list("MAM00001c", "MAR_MEDIUM_EX"))
  met <- data.frame(metabolite_id = rownames(S), compartment = "c", stringsAsFactors = FALSE)
  medium <- data.frame(exchange_reaction_id = "MAR_MEDIUM_EX", lb = -10, ub = 1000, available = TRUE, stringsAsFactors = FALSE)
  gem <- rc_make_gem(S, lb = c(MAR_MEDIUM_EX = -1000), ub = c(MAR_MEDIUM_EX = 1000), metabolite_meta = met)
  gem <- rc_annotate_reaction_roles(gem, medium_table = medium, overwrite_existing = TRUE)
  roles <- stats::setNames(gem$reaction_roles$role, gem$reaction_roles$reaction_id)
  expect_equal(roles[["MAR_MEDIUM_EX"]], "exchange")
})
