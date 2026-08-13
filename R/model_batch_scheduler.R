# Balanced batching across shared structural models.
#
# This overrides the earlier scheduler after microcompass_vmax_cache.R is
# collated. Scheduling changes only task partitioning: every directional target
# remains present exactly once and the numerical LP inputs are unchanged.

.rc_microcompass_vmax_tasks <- function(model_keys, workers) {
  if (is.null(names(model_keys)) || any(!nzchar(names(model_keys)))) {
    stop("Directional vmax model keys require named target rows.",
         call. = FALSE)
  }
  unique_keys <- unique(unname(model_keys))
  n_models <- length(unique_keys)
  workers <- max(1L, as.integer(workers[[1L]]))
  rows_by_model <- lapply(unique_keys, function(model_key) {
    names(model_keys)[model_keys == model_key]
  })
  n_rows <- vapply(rows_by_model, length, integer(1))

  target_batches <- min(
    length(model_keys),
    max(n_models, workers)
  )
  n_batches <- rep.int(1L, n_models)
  remaining <- target_batches - n_models
  while (remaining > 0L) {
    eligible <- which(n_batches < n_rows)
    if (!length(eligible)) break
    load <- n_rows[eligible] / n_batches[eligible]
    chosen <- eligible[[which.max(load)]]
    n_batches[[chosen]] <- n_batches[[chosen]] + 1L
    remaining <- remaining - 1L
  }

  tasks <- list()
  cursor <- 0L
  for (model_index in seq_along(unique_keys)) {
    selected_rows <- rows_by_model[[model_index]]
    batches <- n_batches[[model_index]]
    batch_id <- ceiling(
      seq_along(selected_rows) * batches / length(selected_rows)
    )
    groups <- split(selected_rows, batch_id)
    for (batch_index in seq_along(groups)) {
      cursor <- cursor + 1L
      tasks[[cursor]] <- list(
        model_key = unique_keys[[model_index]],
        row_ids = as.character(groups[[batch_index]])
      )
      names(tasks)[[cursor]] <- paste0(
        "model_", model_index, "__batch_", batch_index
      )
    }
  }
  tasks
}
