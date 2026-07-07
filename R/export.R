#' Export standardized microCOMPASS outputs
#' @export
rc_export_microcompass <- function(result, outdir, write_matrices = TRUE, write_diagnostics = TRUE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (write_matrices) {
    dir.create(file.path(outdir, "04_microcompass"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result$score, file.path(outdir, "04_microcompass", "strict_score_matrix.rds"))
    saveRDS(result$penalty, file.path(outdir, "04_microcompass", "strict_penalty_matrix.rds"))
    saveRDS(result$vmax, file.path(outdir, "04_microcompass", "vmax_matrix.rds"))
    saveRDS(result$feasible, file.path(outdir, "04_microcompass", "feasible_matrix.rds"))
    saveRDS(result$penalty_components, file.path(outdir, "04_microcompass", "penalty_components.rds"))
  }
  if (write_diagnostics) {
    dir.create(file.path(outdir, "02_medium"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "03_microgem"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "04_microcompass"), recursive = TRUE, showWarnings = FALSE)
    if (!is.null(result$medium_scenarios)) utils::write.table(result$medium_scenarios, gzfile(file.path(outdir, "02_medium", "medium_scenarios.tsv.gz"), "wt"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!is.null(result$medium_sensitivity_summary)) utils::write.table(result$medium_sensitivity_summary, gzfile(file.path(outdir, "02_medium", "medium_sensitivity_summary.tsv.gz"), "wt"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!is.null(result$microgem_diagnostics)) utils::write.table(result$microgem_diagnostics, gzfile(file.path(outdir, "03_microgem", "closure_diagnostics.tsv.gz"), "wt"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!is.null(result$microgem_cache_summary)) utils::write.table(result$microgem_cache_summary, gzfile(file.path(outdir, "03_microgem", "microgem_cache_summary.tsv.gz"), "wt"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!is.null(result$lp_diagnostics)) utils::write.table(result$lp_diagnostics, gzfile(file.path(outdir, "04_microcompass", "lp_diagnostics.tsv.gz"), "wt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, "session_info.txt"))
  invisible(outdir)
}
