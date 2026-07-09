#' Write a markdown diagnostics report
#' @export
rc_write_report_md <- function(file,
                               unit_diagnostics = NULL,
                               q95_diagnostics = NULL,
                               gpr_diagnostics = NULL,
                               confidence = NULL,
                               tau_sensitivity = NULL,
                               promiscuity_sensitivity = NULL) {
  lines <- c("# RegCompassR Layer 1 diagnostics", "")
  add_num <- function(label, x) paste0("- ", label, ": ", if (all(is.na(x))) "NA" else signif(stats::median(x, na.rm = TRUE), 4), " median; IQR ", if (all(is.na(x))) "NA" else signif(stats::IQR(x, na.rm = TRUE), 4))
  if (!is.null(unit_diagnostics)) {
    lines <- c(lines, "## Meta-cell diagnostics", paste0("- Meta cells: ", nrow(unit_diagnostics)))
    if ("n_cells" %in% colnames(unit_diagnostics)) lines <- c(lines, add_num("Meta-cell cell count", unit_diagnostics$n_cells))
    if ("low_power_pool" %in% colnames(unit_diagnostics)) lines <- c(lines, paste0("- Low-power meta-cell fraction: ", mean(unit_diagnostics$low_power_pool, na.rm = TRUE)))
    for (nm in intersect(c("pool_seed", "state_source", "state_resolution"), colnames(unit_diagnostics))) lines <- c(lines, paste0("- ", nm, ": ", paste(unique(stats::na.omit(unit_diagnostics[[nm]])), collapse = ", ")))
  }
  if (!is.null(q95_diagnostics)) {
    lines <- c(lines, "", "## Q95 diagnostics", paste0("- Rows: ", nrow(q95_diagnostics)))
    if ("rho_n" %in% colnames(q95_diagnostics)) lines <- c(lines, add_num("Q95 rho_n", q95_diagnostics$rho_n))
    if ("q95_ci_width" %in% colnames(q95_diagnostics)) lines <- c(lines, add_num("Q95 CI width", q95_diagnostics$q95_ci_width))
    if ("q95_power_class" %in% colnames(q95_diagnostics)) {
      tab <- table(q95_diagnostics$q95_power_class, useNA = "ifany")
      lines <- c(lines, paste0("- q95_power_class counts: ", paste(paste(names(tab), as.integer(tab), sep = "="), collapse = ", ")))
    }
    if ("q95_unstable_flag" %in% colnames(q95_diagnostics)) lines <- c(lines, paste0("- q95_unstable_flag fraction: ", mean(q95_diagnostics$q95_unstable_flag, na.rm = TRUE)))
    if ("all_missing_reaction_flag" %in% colnames(q95_diagnostics)) lines <- c(lines, paste0("- all_missing_reaction_flag fraction: ", mean(q95_diagnostics$all_missing_reaction_flag, na.rm = TRUE)))
  }
  if (!is.null(gpr_diagnostics)) {
    lines <- c(lines, "", "## GPR diagnostics", paste0("- Reactions: ", nrow(gpr_diagnostics)))
    if ("missing_subunit_fraction" %in% colnames(gpr_diagnostics)) lines <- c(lines, add_num("Missing subunit fraction", gpr_diagnostics$missing_subunit_fraction))
    if ("missing_subunit_flag" %in% colnames(gpr_diagnostics)) lines <- c(lines, paste0("- Missing-subunit reaction fraction: ", mean(gpr_diagnostics$missing_subunit_flag, na.rm = TRUE)))
  }
  if (!is.null(confidence)) {
    conf <- as.data.frame(confidence)
    lines <- c(lines, "", "## Confidence", paste0("- Rows: ", nrow(conf)))
    if ("reaction_confidence" %in% colnames(conf)) lines <- c(lines, add_num("Reaction confidence", conf$reaction_confidence))
    if ("reaction_confidence_method" %in% colnames(conf)) {
      tab <- table(conf$reaction_confidence_method, useNA = "ifany")
      lines <- c(lines, paste0("- reaction_confidence_method: ", paste(paste(names(tab), as.integer(tab), sep = "="), collapse = ", ")))
    }
    if ("confidence_source" %in% colnames(conf)) {
      tab <- table(conf$confidence_source, useNA = "ifany")
      lines <- c(lines, paste0("- confidence_source counts: ", paste(paste(names(tab), as.integer(tab), sep = "="), collapse = ", ")))
    }
    for (nm in c("no_complete_gpr_group_flag", "reaction_unsupported_by_complete_gpr_flag", "any_incomplete_gpr_group_flag", "low_confidence_reaction_flag")) {
      if (nm %in% colnames(conf)) lines <- c(lines, paste0("- ", nm, " fraction: ", mean(conf[[nm]], na.rm = TRUE)))
    }
    if ("complete_and_group_fraction" %in% colnames(conf)) lines <- c(lines, add_num("Complete AND-group fraction", conf$complete_and_group_fraction))
    if ("best_and_group_observed_fraction" %in% colnames(conf)) lines <- c(lines, add_num("Best AND-group observed fraction", conf$best_and_group_observed_fraction))
  }
  if (!is.null(tau_sensitivity)) lines <- c(lines, "", "## Tau sensitivity", paste0("- Rows: ", nrow(tau_sensitivity)))
  if (!is.null(promiscuity_sensitivity)) lines <- c(lines, "", "## Promiscuity sensitivity", paste0("- Rows: ", nrow(promiscuity_sensitivity)))
  writeLines(lines, con = file)
  invisible(file)
}
