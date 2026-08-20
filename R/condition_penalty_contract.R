# Pando edge eligibility used by the canonical Stage-1 merge.
#
# Conditional GRNs use one frozen exact-edge dictionary. E-star z=0.25 supplies
# continuous condition-specific production coefficients. Formal inference is
# separate and no-fusion; Pando assigns one omnibus P value per exact edge and
# performs BH once across the complete broad-cell-type edge network. RegCompass
# validates and consumes that common topology without recomputing a second gate.
# Target R2 remains diagnostic only on the conditional route.

.RC_PANDO_PENALTY_CORR_THRESHOLD <- 0
.RC_PANDO_PENALTY_ESTIMATE_THRESHOLD <- 0
.RC_PANDO_TARGET_RSQ_THRESHOLD <- 0.05

.rc_target_rsq_threshold <- function(value = getOption(
    "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD)) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) || value < 0 || value >= 1) {
    stop("RegCompass target R2 threshold must be one value in [0, 1).",
         call. = FALSE)
  }
  value[[1L]]
}

.rc_condition_fit_diagnostics_for_coefficients <- function(fit, coefficient) {
  fit_table <- as.data.frame(fit$fit, stringsAsFactors = FALSE)
  required <- c("target", "condition", "fit_status", "rsq")
  if (!all(required %in% colnames(fit_table)) ||
      !all(c("target", "condition") %in% colnames(coefficient))) {
    stop(
      "Condition-GRN diagnostics require target-level fit_status and rsq.",
      call. = FALSE
    )
  }
  fit_key <- paste(
    toupper(trimws(as.character(fit_table$target))),
    as.character(fit_table$condition), sep = "\001"
  )
  coefficient_key <- paste(
    toupper(trimws(as.character(coefficient$target))),
    as.character(coefficient$condition), sep = "\001"
  )
  if (anyNA(fit_key) || anyDuplicated(fit_key)) {
    stop("Condition-GRN target diagnostics are duplicated or incomplete.",
         call. = FALSE)
  }
  index <- match(coefficient_key, fit_key)
  if (anyNA(index)) {
    stop("Condition-GRN coefficients cannot be aligned to diagnostics.",
         call. = FALSE)
  }
  data.frame(
    fit_status = as.character(fit_table$fit_status[index]),
    target_rsq = suppressWarnings(as.numeric(fit_table$rsq[index])),
    stringsAsFactors = FALSE
  )
}

.rc_condition_padj_threshold <- function(fit = NULL, coefficient = NULL) {
  value <- if (!is.null(fit) && !is.null(fit$padj_threshold)) {
    fit$padj_threshold
  } else if (is.data.frame(coefficient) &&
             "padj_threshold" %in% colnames(coefficient)) {
    unique(suppressWarnings(as.numeric(coefficient$padj_threshold)))
  } else {
    0.05
  }
  value <- suppressWarnings(as.numeric(value))
  value <- value[is.finite(value)]
  if (length(value) != 1L || value <= 0 || value >= 1) {
    stop("Condition-GRN padj_threshold must be one value in (0, 1).",
         call. = FALSE)
  }
  value[[1L]]
}

.rc_expected_edge_inference_from_coefficients <- function(coefficient) {
  edge_ids <- unique(as.character(coefficient$edge_id))
  rows <- vector("list", length(edge_ids))
  for (i in seq_along(edge_ids)) {
    edge_id <- edge_ids[[i]]
    index <- which(as.character(coefficient$edge_id) == edge_id)
    valid <- index[
      coefficient$condition_inference_estimable[index] %in% TRUE &
      is.finite(suppressWarnings(as.numeric(
        coefficient$inference_estimate[index]
      ))) &
      is.finite(suppressWarnings(as.numeric(
        coefficient$inference_variance[index]
      ))) &
      suppressWarnings(as.numeric(coefficient$inference_variance[index])) > 0 &
      is.finite(suppressWarnings(as.numeric(coefficient$condition_pval[index])))
    ]
    m <- length(valid)
    statistic <- pval <- NA_real_
    test <- "not_estimable"
    if (m == 1L) {
      one <- valid[[1L]]
      statistic <- suppressWarnings(
        as.numeric(coefficient$inference_statistic[[one]])
      )^2
      pval <- suppressWarnings(as.numeric(coefficient$condition_pval[[one]]))
      test <- "single_condition_exact_t"
    } else if (m > 1L) {
      beta <- suppressWarnings(as.numeric(coefficient$inference_estimate[valid]))
      variance <- suppressWarnings(as.numeric(
        coefficient$inference_variance[valid]
      ))
      statistic <- sum(beta^2 / variance)
      pval <- stats::pchisq(statistic, df = m, lower.tail = FALSE)
      test <- "independent_condition_wald_chisq"
    }
    rows[[i]] <- data.frame(
      edge_id = edge_id,
      edge_df = as.integer(m),
      edge_statistic = statistic,
      edge_pval = pval,
      edge_inference_estimable = m > 0L,
      edge_inference_test = test,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.rc_validate_pando_active_condition_edges <- function(
    coefficient, padj_threshold = NULL) {
  required <- c(
    "edge_id", "target", "condition", "estimate", "penalty_effect",
    "inference_estimate", "inference_se", "inference_variance",
    "inference_statistic", "condition_pval",
    "condition_inference_estimable", "edge_df", "edge_statistic",
    "edge_pval", "edge_padj", "edge_inference_estimable",
    "edge_inference_test", "all_conditions_fit_valid", "edge_supported",
    "statistically_supported", "significant", "pando_estimation_active",
    "active", "active_in_regcompass", "fit_status", "bh_scope",
    "bh_family_size", "penalty_family", "penalty_value", "solver_status",
    "kkt_residual", "iterations"
  )
  if (!is.data.frame(coefficient) ||
      !all(required %in% colnames(coefficient))) {
    stop(
      "Condition-GRN coefficients must retain E-star production, separate ",
      "no-fusion inference, and exact-edge whole-network BH fields.",
      call. = FALSE
    )
  }
  threshold <- if (is.null(padj_threshold)) {
    .rc_condition_padj_threshold(coefficient = coefficient)
  } else {
    value <- suppressWarnings(as.numeric(padj_threshold))
    if (length(value) != 1L || !is.finite(value) ||
        value <= 0 || value >= 1) {
      stop("Condition-GRN padj_threshold must be in (0, 1).",
           call. = FALSE)
    }
    value
  }

  estimate <- suppressWarnings(as.numeric(coefficient$estimate))
  effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  expected_active <- is.finite(estimate)
  if (!identical(
      as.logical(coefficient$pando_estimation_active), expected_active
  ) || !identical(as.logical(coefficient$active), expected_active)) {
    stop("Pando estimation-active flags must track finite E-star coefficients.",
         call. = FALSE)
  }
  comparable <- is.finite(estimate) & is.finite(effect)
  if (any(is.finite(estimate) != is.finite(effect)) ||
      any(abs(estimate[comparable] - effect[comparable]) > 1e-12)) {
    stop("penalty_effect must equal the production E-star coefficient.",
         call. = FALSE)
  }

  if (any(as.character(coefficient$bh_scope) !=
          "exact_edge_whole_cell_type_network_BH")) {
    stop("Conditional topology must use exact-edge whole-network BH.",
         call. = FALSE)
  }

  expected <- .rc_expected_edge_inference_from_coefficients(coefficient)
  observed_edge <- coefficient[!duplicated(as.character(coefficient$edge_id)),
                               , drop = FALSE]
  observed_edge <- observed_edge[
    match(expected$edge_id, as.character(observed_edge$edge_id)), , drop = FALSE
  ]
  if (anyNA(observed_edge$edge_id)) {
    stop("Exact-edge inference cannot be aligned to coefficient rows.",
         call. = FALSE)
  }
  compare_numeric <- function(observed, wanted, tol = 1e-10) {
    observed <- suppressWarnings(as.numeric(observed))
    wanted <- suppressWarnings(as.numeric(wanted))
    finite <- is.finite(observed) & is.finite(wanted)
    !any(is.finite(observed) != is.finite(wanted)) &&
      !any(abs(observed[finite] - wanted[finite]) > tol)
  }
  if (!identical(as.integer(observed_edge$edge_df), expected$edge_df) ||
      !identical(as.logical(observed_edge$edge_inference_estimable),
                 expected$edge_inference_estimable) ||
      !identical(as.character(observed_edge$edge_inference_test),
                 expected$edge_inference_test) ||
      !compare_numeric(observed_edge$edge_statistic,
                       expected$edge_statistic) ||
      !compare_numeric(observed_edge$edge_pval, expected$edge_pval)) {
    stop("Stored exact-edge omnibus inference is inconsistent.",
         call. = FALSE)
  }

  expected_padj <- rep(NA_real_, nrow(expected))
  valid_edge <- which(
    expected$edge_inference_estimable %in% TRUE &
      is.finite(expected$edge_pval)
  )
  if (length(valid_edge)) {
    expected_padj[valid_edge] <- stats::p.adjust(
      expected$edge_pval[valid_edge], method = "BH"
    )
  }
  if (!compare_numeric(observed_edge$edge_padj, expected_padj)) {
    stop("Stored edge_padj is not BH across the complete exact-edge network.",
         call. = FALSE)
  }
  if (any(as.integer(observed_edge$bh_family_size) != length(valid_edge))) {
    stop("Stored exact-edge BH family size is inconsistent.", call. = FALSE)
  }

  edge_ids <- as.character(expected$edge_id)
  for (edge_id in edge_ids) {
    index <- which(as.character(coefficient$edge_id) == edge_id)
    valid_production <- all(
      as.character(coefficient$fit_status[index]) == "ok" &
        is.finite(effect[index])
    )
    unique_padj <- unique(suppressWarnings(as.numeric(
      coefficient$edge_padj[index]
    )))
    unique_support <- unique(as.logical(coefficient$edge_supported[index]))
    unique_active <- unique(as.logical(
      coefficient$active_in_regcompass[index]
    ))
    unique_valid <- unique(as.logical(
      coefficient$all_conditions_fit_valid[index]
    ))
    if (length(unique_padj) != 1L || length(unique_support) != 1L ||
        length(unique_active) != 1L || length(unique_valid) != 1L) {
      stop("Exact-edge topology fields must be identical across conditions.",
           call. = FALSE)
    }
    expected_support <- valid_production && is.finite(unique_padj) &&
      unique_padj < threshold
    if (!identical(unique_valid[[1L]], valid_production) ||
        !identical(unique_support[[1L]], expected_support) ||
        !identical(unique_active[[1L]], expected_support) ||
        any(as.logical(coefficient$statistically_supported[index]) !=
            expected_support) ||
        any(as.logical(coefficient$significant[index]) != expected_support)) {
      stop("Pando common exact-edge topology flags are inconsistent.",
           call. = FALSE)
    }
    generic_p <- suppressWarnings(as.numeric(coefficient$pval[index]))
    generic_q <- suppressWarnings(as.numeric(coefficient$padj[index]))
    if (any(is.finite(generic_p) !=
            is.finite(as.numeric(coefficient$edge_pval[index]))) ||
        any(is.finite(generic_q) !=
            is.finite(as.numeric(coefficient$edge_padj[index]))) ||
        any(abs(generic_p[is.finite(generic_p)] -
                as.numeric(coefficient$edge_pval[index][is.finite(generic_p)])) >
            1e-12) ||
        any(abs(generic_q[is.finite(generic_q)] -
                as.numeric(coefficient$edge_padj[index][is.finite(generic_q)])) >
            1e-12)) {
      stop("Generic pval/padj fields must mirror exact-edge inference.",
           call. = FALSE)
    }
  }

  family <- as.character(coefficient$penalty_family)
  value <- suppressWarnings(as.numeric(coefficient$penalty_value))
  solver <- as.character(coefficient$solver_status)
  kkt <- suppressWarnings(as.numeric(coefficient$kkt_residual))
  iteration <- suppressWarnings(as.integer(coefficient$iterations))
  if (anyNA(family) ||
      any(family != "information_scaled_sparse_deviation") ||
      any(!is.finite(value)) || any(abs(value - 0.25) > 1e-15) ||
      anyNA(solver) || any(solver != "ok") ||
      any(!is.finite(kkt)) || any(kkt < 0) ||
      anyNA(iteration) || any(iteration < 0L)) {
    stop("Conditional coefficients are not converged E-star z=0.25 results.",
         call. = FALSE)
  }
  invisible(as.logical(coefficient$active_in_regcompass))
}

.rc_condition_target_rsq <- function(coefficient) {
  if ("target_rsq" %in% colnames(coefficient)) {
    return(suppressWarnings(as.numeric(coefficient$target_rsq)))
  }
  if ("rsq" %in% colnames(coefficient)) {
    return(suppressWarnings(as.numeric(coefficient$rsq)))
  }
  stop("Condition-GRN diagnostics require full-data target R2 (`rsq`).",
       call. = FALSE)
}

.rc_condition_penalty_gate <- function(
    coefficient, padj_threshold = NULL, target_rsq_threshold = NULL) {
  threshold <- if (is.null(padj_threshold)) {
    .rc_condition_padj_threshold(coefficient = coefficient)
  } else {
    .rc_condition_padj_threshold(
      coefficient = transform(coefficient, padj_threshold = padj_threshold)
    )
  }
  invisible(.rc_target_rsq_threshold(
    target_rsq_threshold %||% getOption(
      "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD
    )
  ))
  gate <- .rc_validate_pando_active_condition_edges(
    coefficient, padj_threshold = threshold
  )
  effect <- suppressWarnings(as.numeric(coefficient$penalty_effect))
  gate & as.character(coefficient$fit_status) == "ok" & is.finite(effect)
}

.rc_apply_condition_penalty_gate <- function(
    fit, target_rsq_threshold = NULL) {
  .rc_require_pando_condition_grn_fit(fit)
  threshold <- .rc_condition_padj_threshold(fit = fit)
  rsq_threshold <- .rc_target_rsq_threshold(
    target_rsq_threshold %||% getOption(
      "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD
    )
  )
  coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
  diagnostics <- .rc_condition_fit_diagnostics_for_coefficients(
    fit, coefficient
  )
  coefficient$fit_status <- diagnostics$fit_status
  coefficient$target_rsq <- diagnostics$target_rsq
  coefficient$rsq <- diagnostics$target_rsq
  coefficient$padj_threshold <- threshold
  coefficient$target_rsq_threshold <- rsq_threshold
  coefficient$target_model_supported <-
    coefficient$fit_status == "ok" & is.finite(coefficient$target_rsq) &
    coefficient$target_rsq >= rsq_threshold
  gate <- .rc_condition_penalty_gate(
    coefficient,
    padj_threshold = threshold,
    target_rsq_threshold = rsq_threshold
  )
  coefficient$penalty_eligible <- gate
  coefficient$active_in_condition <- gate
  fit$coefficients <- coefficient
  fit$regcompass_penalty_filter <- paste0(
    "Pando common exact edge: whole-network BH edge_padj < ",
    format(threshold, trim = TRUE),
    " with valid finite production beta_E in every condition"
  )
  fit$regcompass_fit_status_filter <- "fit_status == 'ok' in every condition"
  fit$regcompass_target_rsq_filter <- paste0(
    "diagnostic only: rsq >= ", format(rsq_threshold, trim = TRUE)
  )
  fit$regcompass_target_rsq_definition <-
    "scheme_e_z025_full_data_R2_diagnostic"
  fit$regcompass_significance_role <-
    "consume_Pando_exact_edge_whole_network_BH_common_topology;R2_diagnostic_only"
  fit$regcompass_padj_threshold <- threshold
  fit$regcompass_target_rsq_threshold <- rsq_threshold
  fit
}

.rc_filter_standard_pando_edges <- function(
    table, padj_threshold = .rc_standard_pando_padj_default,
    target_rsq_threshold = NULL) {
  required <- c("estimate", "padj", "rsq")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop(
      "Standard Pando requires estimate, padj and full-data rsq columns for ",
      "RegCompass penalty filtering.", call. = FALSE
    )
  }
  threshold <- suppressWarnings(as.numeric(padj_threshold))
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold >= 1) {
    stop("Standard Pando `padj_threshold` must be one value in (0, 1).",
         call. = FALSE)
  }
  rsq_threshold <- .rc_target_rsq_threshold(
    target_rsq_threshold %||% getOption(
      "RegCompassR.target_rsq_threshold", .RC_PANDO_TARGET_RSQ_THRESHOLD
    )
  )
  estimate <- suppressWarnings(as.numeric(table$estimate))
  padj <- suppressWarnings(as.numeric(table$padj))
  rsq <- suppressWarnings(as.numeric(table$rsq))
  keep <- is.finite(estimate) & is.finite(padj) & padj < threshold &
    is.finite(rsq) & rsq >= rsq_threshold
  if ("estimable" %in% colnames(table)) {
    keep <- keep & table$estimable %in% TRUE
  }
  answer <- table[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    estimable = if ("estimable" %in% colnames(table)) TRUE else NA,
    padj = paste0("< ", format(threshold, trim = TRUE)),
    rsq = paste0(">= ", format(rsq_threshold, trim = TRUE)),
    rsq_definition = "selected_lambda_full_data_R2"
  )
  answer
}
