#' Export standardized microCOMPASS outputs
#' @export
rc_export_microcompass <- function(result, outdir,
                                   write_matrices = TRUE,
                                   write_diagnostics = TRUE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  microcompass_dir <- file.path(outdir, "04_microcompass")
  model_dir <- file.path(outdir, "03_models")

  if (isTRUE(write_matrices)) {
    dir.create(
      microcompass_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    saveRDS(
      result$score,
      file.path(microcompass_dir, "strict_score_matrix.rds")
    )
    saveRDS(
      result$penalty,
      file.path(microcompass_dir, "strict_penalty_matrix.rds")
    )
    saveRDS(
      result$vmax,
      file.path(microcompass_dir, "vmax_matrix.rds")
    )
    saveRDS(
      result$feasible,
      file.path(microcompass_dir, "feasible_matrix.rds")
    )
    if (!is.null(result$evaluated)) {
      saveRDS(
        result$evaluated,
        file.path(microcompass_dir, "evaluated_matrix.rds")
      )
    }
    saveRDS(
      result$penalty_components,
      file.path(microcompass_dir, "penalty_components.rds")
    )
    saveRDS(
      result$params,
      file.path(microcompass_dir, "run_parameters.rds")
    )
  }

  if (isTRUE(write_diagnostics)) {
    medium_dir <- file.path(outdir, "02_medium")
    dir.create(medium_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(
      microcompass_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (!is.null(result$medium_scenarios)) {
      utils::write.table(
        result$medium_scenarios,
        gzfile(
          file.path(medium_dir, "medium_scenarios.tsv.gz"),
          "wt"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
    if (!is.null(result$model_cache_summary)) {
      utils::write.table(
        result$model_cache_summary,
        gzfile(
          file.path(model_dir, "model_cache_summary.tsv.gz"),
          "wt"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
    if (!is.null(result$model_diagnostics)) {
      utils::write.table(
        result$model_diagnostics,
        gzfile(
          file.path(model_dir, "model_diagnostics.tsv.gz"),
          "wt"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
    if (!is.null(result$lp_diagnostics)) {
      utils::write.table(
        result$lp_diagnostics,
        gzfile(
          file.path(microcompass_dir, "lp_diagnostics.tsv.gz"),
          "wt"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(outdir, "session_info.txt")
  )
  invisible(outdir)
}
