#' Write a minimal markdown diagnostics report
#' @export
rc_write_report_md <- function(file,
                               pool_diagnostics = NULL,
                               q95_diagnostics = NULL,
                               gpr_diagnostics = NULL,
                               confidence = NULL) {
  lines <- c("# RegCompassR Layer 1 diagnostics", "")
  if (!is.null(pool_diagnostics)) {
    lines <- c(lines, "## Pool diagnostics", paste0("- Pools: ", nrow(pool_diagnostics)))
    if ("low_power_pool" %in% colnames(pool_diagnostics)) lines <- c(lines, paste0("- Low-power pool fraction: ", mean(pool_diagnostics$low_power_pool, na.rm = TRUE)))
  }
  if (!is.null(q95_diagnostics)) {
    lines <- c(lines, "", "## Q95 diagnostics", paste0("- Rows: ", nrow(q95_diagnostics)))
    for (nm in intersect(c("q95_low_power", "q95_very_low_power"), colnames(q95_diagnostics))) lines <- c(lines, paste0("- ", nm, " fraction: ", mean(q95_diagnostics[[nm]], na.rm = TRUE)))
  }
  if (!is.null(gpr_diagnostics)) {
    lines <- c(lines, "", "## GPR diagnostics", paste0("- Reactions: ", nrow(gpr_diagnostics)))
    if ("missing_subunit_flag" %in% colnames(gpr_diagnostics)) lines <- c(lines, paste0("- Missing-subunit reaction fraction: ", mean(gpr_diagnostics$missing_subunit_flag, na.rm = TRUE)))
  }
  if (!is.null(confidence)) {
    lines <- c(lines, "", "## Confidence", paste0("- Rows: ", nrow(as.data.frame(confidence))))
  }
  writeLines(lines, con = file)
  invisible(file)
}
