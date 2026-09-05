# Exact directional COMPASS Step 2 compiler.
#
# This implementation corrects the zero-length directional-ID edge case in
# the PR #338 compiler. Base R paste0() recycles a non-empty suffix against a
# zero-length character vector unless recycle0 = TRUE, which can create a
# spurious "::forward" or "::reverse" name when an entire model has no flux in
# that direction. The explicit length guards below preserve the exact number
# of directional variables and do not alter the LP feasible set.
#
# The uniquely named implementation is bound to the production symbol at the
# end of this file. This keeps the API surface auditable: there is only one
# direct top-level function definition for each implementation name.

.rc_compass_step2_prepare_exact <- function(
    S, lb, ub, target_reaction, vmax_result,
    target_direction = c("forward", "reverse"),
    omega = 0.95, flux_threshold = 1e-8) {
  target_direction <- match.arg(target_direction)
  if (!is.numeric(omega) || length(omega) != 1L ||
      !is.finite(omega) || omega <= 0 || omega > 1) {
    stop("`omega` must be in (0, 1].", call. = FALSE)
  }
  if (!is.list(vmax_result) ||
      !all(c("feasible", "vmax", "status") %in% names(vmax_result))) {
    stop("`vmax_result` must be a directional Step 1 result.", call. = FALSE)
  }
  reactions <- colnames(S)
  if (is.null(reactions) || anyNA(reactions) || any(!nzchar(reactions)) ||
      anyDuplicated(reactions)) {
    stop("`S` must have unique non-empty reaction IDs in colnames().",
         call. = FALSE)
  }
  if (!target_reaction %in% reactions) {
    stop("`target_reaction` is missing from the stoichiometric matrix.",
         call. = FALSE)
  }
  lb <- rc_align_bound(lb, reactions, default = -1000, name = "lb")
  ub <- rc_align_bound(ub, reactions, default = 1000, name = "ub")
  if (any(lb > ub)) {
    stop("Reaction lower bounds cannot exceed upper bounds.", call. = FALSE)
  }
  early <- if (!isTRUE(vmax_result$feasible)) {
    list(
      feasible = FALSE,
      penalty = NA_real_,
      vmax = as.numeric(vmax_result$vmax),
      solver_status = as.character(vmax_result$status),
      step1_status = as.character(vmax_result$status),
      step2_status = "not_run",
      solver_backend = "not_run",
      flux = numeric()
    )
  } else {
    NULL
  }
  if (!is.null(early)) {
    return(list(runnable = FALSE, reactions = reactions, result = early))
  }

  vmax <- as.numeric(vmax_result$vmax)
  if (length(vmax) != 1L || !is.finite(vmax) || vmax < flux_threshold) {
    stop("Cached directional vmax is not a positive feasible value.",
         call. = FALSE)
  }

  S <- .rc_as_dgCMatrix(S)
  n_reactions <- ncol(S)
  target_index <- match(target_reaction, reactions)

  # COMPASS Step 2 uses exact non-negative directional variables. CORDA2 has
  # its own structural directional split and directional Vmax remains signed.
  forward_index <- which(ub > 0)
  reverse_index <- which(lb < 0)
  direction_reaction_index <- c(forward_index, reverse_index)
  direction_sign <- c(
    rep.int(1, length(forward_index)),
    rep.int(-1, length(reverse_index))
  )
  if (!length(direction_reaction_index)) {
    stop("The completed GEM contains no flux-carrying reaction direction.",
         call. = FALSE)
  }

  forward_id <- if (length(forward_index)) {
    paste0(reactions[forward_index], "::forward")
  } else {
    character()
  }
  reverse_id <- if (length(reverse_index)) {
    paste0(reactions[reverse_index], "::reverse")
  } else {
    character()
  }
  variable_id <- c(forward_id, reverse_id)
  if (length(variable_id) != length(direction_reaction_index)) {
    stop("Directional Step 2 variable identifiers are malformed.",
         call. = FALSE)
  }

  forward_S <- S[, forward_index, drop = FALSE]
  reverse_S <- S[, reverse_index, drop = FALSE]
  if (length(reverse_S@x)) reverse_S@x <- -reverse_S@x
  A <- cbind(forward_S, reverse_S)
  if (ncol(A) != length(variable_id)) {
    stop("Directional Step 2 matrix and variable IDs do not align.",
         call. = FALSE)
  }
  colnames(A) <- variable_id
  rm(forward_S, reverse_S)

  lower <- c(
    pmax(lb[forward_index], 0),
    pmax(-ub[reverse_index], 0)
  )
  upper <- c(
    ub[forward_index],
    -lb[reverse_index]
  )
  names(lower) <- names(upper) <- variable_id

  target_sign <- if (identical(target_direction, "forward")) 1 else -1
  target_variable <- which(
    direction_reaction_index == target_index & direction_sign == target_sign
  )
  if (length(target_variable) != 1L) {
    stop(
      "Cached directional vmax is feasible but the requested target direction ",
      "is absent from the exact Step 2 directional model.",
      call. = FALSE
    )
  }
  opposite_variable <- which(
    direction_reaction_index == target_index & direction_sign == -target_sign
  )
  required_flux <- omega * vmax
  bound_tolerance <- max(
    flux_threshold,
    32 * .Machine$double.eps * max(
      1, abs(required_flux), abs(upper[[target_variable]])
    )
  )
  if (required_flux > upper[[target_variable]] + bound_tolerance) {
    stop(
      "Cached directional vmax is inconsistent with the completed-GEM target ",
      "bound after exact directional compilation.",
      call. = FALSE
    )
  }
  lower[[target_variable]] <- max(
    lower[[target_variable]],
    min(required_flux, upper[[target_variable]])
  )
  if (length(opposite_variable)) {
    if (lower[[opposite_variable]] > bound_tolerance) {
      stop(
        "The opposite target direction has a positive compulsory lower bound; ",
        "this contradicts the feasible signed directional Vmax contract.",
        call. = FALSE
      )
    }
    lower[[opposite_variable]] <- 0
    upper[[opposite_variable]] <- 0
  }

  list(
    runnable = TRUE,
    reactions = reactions,
    template = list(
      A = A,
      lhs = rep(0, nrow(A)),
      rhs = rep(0, nrow(A)),
      lb = as.numeric(lower),
      ub = as.numeric(upper),
      n_reactions = n_reactions,
      n_variables = length(direction_reaction_index),
      reactions = reactions,
      variable_id = variable_id,
      direction_reaction_index = as.integer(direction_reaction_index),
      direction_sign = as.numeric(direction_sign),
      target_variable_index = as.integer(target_variable),
      opposite_variable_index = as.integer(opposite_variable),
      required_flux = as.numeric(required_flux),
      formulation = "compass_directional_nonnegative_exact_v1",
      vmax = vmax,
      step1_status = as.character(vmax_result$status),
      target_reaction = target_reaction,
      target_direction = target_direction
    )
  )
}

# Production binding. This is intentionally a plain alias, not a saved
# legacy/base/implementation copy: callers continue to use the established
# private symbol while the exact corrected compiler has one auditable definition.
.rc_compass_step2_prepare <- .rc_compass_step2_prepare_exact
