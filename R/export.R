#' Export standardized microCOMPASS outputs
#' @export
rc_export_microcompass <- function(result, outdir,
                                   write_matrices = TRUE,
                                   write_diagnostics = TRUE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  microcompass_dir <- file.path(outdir, "04_microcompass")
  model_dir <- file.path(outdir, "03_models")

  write_tsv_gz <- function(value, path) {
    connection <- gzfile(path, "wt")
    on.exit(close(connection), add = TRUE)
    utils::write.table(
      value,
      connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    invisible(path)
  }

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
      write_tsv_gz(
        result$medium_scenarios,
        file.path(medium_dir, "medium_scenarios.tsv.gz")
      )
    }
    if (!is.null(result$model_cache_summary)) {
      write_tsv_gz(
        result$model_cache_summary,
        file.path(model_dir, "model_cache_summary.tsv.gz")
      )
    }
    if (!is.null(result$model_diagnostics)) {
      write_tsv_gz(
        result$model_diagnostics,
        file.path(model_dir, "model_diagnostics.tsv.gz")
      )
    }
    if (!is.null(result$lp_diagnostics)) {
      write_tsv_gz(
        result$lp_diagnostics,
        file.path(microcompass_dir, "lp_diagnostics.tsv.gz")
      )
    }
  }

  writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(outdir, "session_info.txt")
  )
  invisible(outdir)
}
