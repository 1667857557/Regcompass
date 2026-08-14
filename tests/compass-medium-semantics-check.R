suppressPackageStartupMessages(library(Matrix))

`%||%` <- function(x, y) if (is.null(x)) y else x

rc_validate_gem <- function(gem) {
  S <- methods::as(gem$S, "dgCMatrix")
  reactions <- colnames(S)
  metabolites <- rownames(S)
  lb <- as.numeric(gem$lb[reactions])
  ub <- as.numeric(gem$ub[reactions])
  names(lb) <- reactions
  names(ub) <- reactions
  stopifnot(!anyNA(lb), !anyNA(ub), !any(lb > ub))
  list(
    S = S, lb = lb, ub = ub,
    reactions = reactions, metabolites = metabolites
  )
}

rc_annotate_reaction_roles <- function(gem, medium_table = NULL) gem

source("R/full_gem.R")
source("R/medium.R")
source("R/compass_medium_semantics.R")
source("R/layer2_corda_parent_contract.R")

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
gem <- list(
  S = S,
  lb = stats::setNames(rep(-1000, length(reactions)), reactions),
  ub = stats::setNames(rep(1000, length(reactions)), reactions),
  reaction_meta = data.frame(
    reaction_id = reactions,
    role = c(rep("exchange", 5), "internal"),
    stringsAsFactors = FALSE
  )
)

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
stopifnot(
  identical(unname(applied$gem$lb[["EX_reverse"]]), -10),
  identical(unname(applied$gem$ub[["EX_reverse"]]), 1000),
  identical(unname(applied$gem$lb[["EX_forward"]]), -1000),
  identical(unname(applied$gem$ub[["EX_forward"]]), 4),
  identical(unname(applied$gem$lb[["EX_unlisted"]]), -1),
  identical(unname(applied$gem$ub[["EX_unlisted"]]), 1000),
  identical(unname(applied$gem$lb[["EX_close_reverse"]]), 0),
  identical(unname(applied$gem$ub[["EX_close_reverse"]]), 1000),
  identical(unname(applied$gem$lb[["EX_close_forward"]]), -1000),
  identical(unname(applied$gem$ub[["EX_close_forward"]]), 0),
  identical(unname(applied$gem$lb[["R_internal"]]), -1000),
  identical(unname(applied$gem$ub[["R_internal"]]), 1000),
  identical(
    applied$gem$medium_semantics_version,
    "compass_exchange_bounds_v2"
  )
)

contract <- .rc_assert_medium_bounds_only(gem, applied$gem, medium)
stopifnot(
  "EX_unlisted" %in% contract$changed_unlisted_exchange_reactions,
  contract$n_removed_reactions == 0L
)

invalid <- applied$gem
invalid$ub[["R_internal"]] <- 5
invalid_contract <- try(
  .rc_assert_medium_bounds_only(gem, invalid, medium), silent = TRUE
)
stopifnot(inherits(invalid_contract, "try-error"))

custom <- data.frame(
  medium_scenario_id = "custom",
  exchange_reaction_id = "EX_reverse",
  condition = "all",
  lb = -0.2,
  ub = 5,
  available = TRUE,
  stringsAsFactors = FALSE
)
custom_applied <- rc_apply_medium_constraints(gem, custom)
stopifnot(
  identical(unname(custom_applied$gem$lb[["EX_reverse"]]), -0.2),
  identical(unname(custom_applied$gem$ub[["EX_reverse"]]), 5)
)

model_bounds <- .rc_compass_model_bound_medium(gem, exchange_limit = 1)
reverse_row <- model_bounds[
  model_bounds$exchange_reaction_id == "EX_reverse", , drop = FALSE
]
forward_row <- model_bounds[
  model_bounds$exchange_reaction_id == "EX_forward", , drop = FALSE
]
stopifnot(
  reverse_row$lb == -1,
  reverse_row$ub == 1000,
  forward_row$lb == -1000,
  forward_row$ub == 1,
  all(model_bounds$bound_scope == "uptake")
)

full <- rc_build_full_gem(gem, medium_table = medium)
corda_parent <- .rc_corda_parent(
  gem,
  medium_table = medium,
  forbidden_roles = character(),
  solver = "highs",
  time_limit = 60
)
stopifnot(
  identical(colnames(full$S), reactions),
  identical(colnames(corda_parent$S), reactions),
  identical(unname(full$lb[["EX_unlisted"]]), -1),
  identical(unname(corda_parent$lb[["EX_unlisted"]]), -1),
  identical(unname(full$ub[["EX_reverse"]]), 1000),
  identical(unname(corda_parent$ub[["EX_forward"]]), 4),
  identical(full$build_params$n_medium_removed_reactions, 0L),
  identical(corda_parent$corda_parent_prepruning, "none"),
  identical(corda_parent$corda_parent_role_blocking, "none")
)

identity_1 <- .rc_full_gem_resolved_bounds_identity(gem, full)
medium_provenance_only <- medium
medium_provenance_only$concentration_mM <- c(99, 88, 77, 66)
full_same_bounds <- rc_build_full_gem(
  gem, medium_table = medium_provenance_only
)
identity_same <- .rc_full_gem_resolved_bounds_identity(gem, full_same_bounds)
medium_changed_bounds <- medium
medium_changed_bounds$exchange_limit[[1L]] <- 11
full_changed_bounds <- rc_build_full_gem(
  gem, medium_table = medium_changed_bounds
)
identity_changed <- .rc_full_gem_resolved_bounds_identity(
  gem, full_changed_bounds
)
stopifnot(
  is.character(identity_1$fingerprint), nzchar(identity_1$fingerprint),
  identical(identity_1$fingerprint, identity_same$fingerprint),
  !identical(identity_1$fingerprint, identity_changed$fingerprint),
  identity_1$n_changed_bounds_vs_reference > 0L,
  !identity_1$resolved_bounds_identical_to_reference,
  identical(
    .rc_full_gem_medium_fingerprint(gem, full),
    identity_1$fingerprint
  )
)

message("COMPASS exchange-medium and Layer 2 handoff regression passed.")
