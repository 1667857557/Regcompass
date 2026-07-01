#' Write a minimal Layer 1 diagnostics report
#' @export
rc_write_report_md <- function(file,
                               layer1 = NULL,
                               pool_diagnostics = NULL,
                               q95_diagnostics = NULL,
                               gpr_diagnostics = NULL,
                               reaction_confidence = NULL,
                               tau_sensitivity = NULL) {
  if (!is.null(layer1) && is.list(layer1$minimal_diagnostics)) {
    pool_diagnostics <- layer1$minimal_diagnostics$pool_diagnostics
    q95_diagnostics <- layer1$minimal_diagnostics$q95_diagnostics
    gpr_diagnostics <- layer1$minimal_diagnostics$gpr_diagnostics
    tau_sensitivity <- layer1$minimal_diagnostics$tau_sensitivity
    reaction_confidence <- layer1$reaction_confidence
  }

  add_num <- function(label, x) paste0("- ", label, ": ", if (all(is.na(x))) "NA" else signif(stats::median(x, na.rm = TRUE), 4),
                                      " median; IQR ", if (all(is.na(x))) "NA" else signif(stats::IQR(x, na.rm = TRUE), 4))
  lines <- c("# RegCompassR simplified Layer 1 report", "",
             "Main model: RNA raw counts -> pool sum -> log2(CPM + 1) -> robust z -> sigmoid -> sqrt promiscuity -> Boltzmann AND tau=0.20 -> OR sum -> cell-type Q95 shrinkage.",
             "")

  if (!is.null(pool_diagnostics)) {
    lines <- c(lines, "## Pool diagnostics", paste0("- Pools: ", nrow(pool_diagnostics)))
    if ("n_cells" %in% colnames(pool_diagnostics)) lines <- c(lines, add_num("Pool cell count", pool_diagnostics$n_cells))
    if ("low_power_pool" %in% colnames(pool_diagnostics)) lines <- c(lines, paste0("- Low-power pool fraction: ", signif(mean(pool_diagnostics$low_power_pool, na.rm = TRUE), 4)))
    if ("GPR_gene_detection_rate" %in% colnames(pool_diagnostics)) lines <- c(lines, add_num("GPR gene detection rate", pool_diagnostics$GPR_gene_detection_rate))
  }
  if (!is.null(q95_diagnostics)) {
    lines <- c(lines, "", "## Q95 diagnostics", paste0("- Rows: ", nrow(q95_diagnostics)))
    if ("rho_n" %in% colnames(q95_diagnostics)) lines <- c(lines, add_num("Q95 shrinkage rho", q95_diagnostics$rho_n))
    if ("q95_low_power" %in% colnames(q95_diagnostics)) lines <- c(lines, paste0("- q95_low_power fraction: ", signif(mean(q95_diagnostics$q95_low_power, na.rm = TRUE), 4)))
  }
  if (!is.null(gpr_diagnostics)) {
    lines <- c(lines, "", "## GPR diagnostics", paste0("- Reactions: ", nrow(gpr_diagnostics)))
    if ("missing_gpr_gene_fraction" %in% colnames(gpr_diagnostics)) lines <- c(lines, add_num("Missing GPR gene fraction", gpr_diagnostics$missing_gpr_gene_fraction))
  }
  if (!is.null(reaction_confidence)) {
    lines <- c(lines, "", "## Reaction confidence", paste0("- Rows: ", nrow(reaction_confidence)))
    if ("reaction_confidence" %in% colnames(reaction_confidence)) lines <- c(lines, add_num("Reaction confidence", reaction_confidence$reaction_confidence))
  }
  if (!is.null(tau_sensitivity)) {
    lines <- c(lines, "", "## Tau sensitivity")
    if ("tau_sensitive_flag" %in% colnames(tau_sensitivity)) lines <- c(lines, paste0("- tau_sensitive reaction fraction: ", signif(mean(tau_sensitivity$tau_sensitive_flag, na.rm = TRUE), 4)))
  }

  writeLines(lines, con = file)
  invisible(file)
}
