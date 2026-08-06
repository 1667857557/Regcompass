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

.rc_medium_exchange_metabolites <- function(gem, exchange_meta, validated) {
  data.frame(
    exchange_reaction_id = as.character(exchange_meta$reaction_id),
    metabolite_id = paste0("m_", seq_len(nrow(exchange_meta))),
    gem_metabolite_name = as.character(exchange_meta$reaction_id),
    mapping_source = "synthetic",
    stringsAsFactors = FALSE
  )
}

source("R/compass_medium_semantics.R")

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

fingerprint_1 <- .rc_full_gem_medium_fingerprint(medium)
medium_2 <- medium
medium_2$exchange_limit[[1L]] <- 11
fingerprint_2 <- .rc_full_gem_medium_fingerprint(medium_2)
stopifnot(
  is.character(fingerprint_1), nzchar(fingerprint_1),
  !identical(fingerprint_1, fingerprint_2)
)

message("COMPASS exchange-medium semantics regression passed.")
