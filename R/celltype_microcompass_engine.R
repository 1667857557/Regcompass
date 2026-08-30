# Directional LP scoring on cell-type-specific structural models.

.rc_celltype_model_contract <- function(model_cache) {
  if (!is.list(model_cache) || !length(model_cache)) {
    stop("The cell-type model cache is empty.", call. = FALSE)
  }
  files <- vapply(model_cache, function(entry) as.character(entry$file), character(1))
  representative <- match(unique(files), files)
  rows <- lapply(representative, function(i) {
    entry <- model_cache[[i]]
    model <- .rc_read_celltype_union_gem(
      entry$file, entry$cell_type, entry$medium_scenario,
      entry$file_checksum
    )
    validated <- rc_validate_gem(model)
    data.frame(
      cell_type = as.character(entry$cell_type),
      medium_scenario = as.character(entry$medium_scenario),
      model_file = as.character(entry$file),
      model_file_checksum = unname(tools::md5sum(entry$file)[[1L]]),
      n_reactions = ncol(validated$S),
      n_metabolites = nrow(validated$S),
      reaction_order_checksum =
        .rc_microcompass_object_checksum(colnames(validated$S)),
      metabolite_order_checksum =
        .rc_microcompass_object_checksum(rownames(validated$S)),
      stoichiometry_bounds_checksum =
        .rc_microcompass_object_checksum(list(
          S = validated$S, lb = validated$lb, ub = validated$ub
        )),
      shared_across_conditions = TRUE,
      shared_across_cell_types = FALSE,
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer[order(answer$cell_type, answer$medium_scenario), , drop = FALSE]
}
