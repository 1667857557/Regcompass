.rc_bind_frames_fill <- function(x) {
  x <- x[vapply(
    x,
    function(z) is.data.frame(z) && nrow(z) > 0L,
    logical(1)
  )]
  if (!length(x)) return(data.frame())
  columns <- unique(unlist(lapply(x, colnames), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(columns, colnames(z))
    for (column in missing) z[[column]] <- NA
    z[, columns, drop = FALSE]
  })
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.rc_metacell_logcpm <- function(
    counts, scale_factor = 1e6, library_size = NULL) {
  counts <- .rc_as_dgCMatrix(counts)
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("`scale_factor` must be one positive finite number.", call. = FALSE)
  }
  normalization_scope <- "input_matrix_library_size"
  if (is.null(library_size)) {
    library_size <- Matrix::colSums(counts)
  } else {
    normalization_scope <- "full_transcriptome_library_size_before_gpr_filter"
    if (!is.null(names(library_size))) {
      missing <- setdiff(colnames(counts), names(library_size))
      if (length(missing)) {
        stop(
          "`library_size` is missing metacells: ",
          paste(utils::head(missing, 10L), collapse = ", "),
          call. = FALSE
        )
      }
      library_size <- library_size[colnames(counts)]
    }
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) ||
      any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop(
      "`library_size` must contain one positive finite value per metacell.",
      call. = FALSE
    )
  }
  scaled <- counts %*% Matrix::Diagonal(x = scale_factor / library_size)
  answer <- log1p(scaled)
  dimnames(answer) <- dimnames(counts)
  attr(answer, "normalization_scope") <- normalization_scope
  attr(answer, "library_size") <- stats::setNames(
    library_size,
    colnames(counts)
  )
  answer
}

.rc_merge_stratum_meta_modules <- function(artifacts) {
  names_to_merge <- c(
    "sample_status", "tf_peak_gene_all", "tf_peak_gene_significant",
    "metabolic_gene_nodes", "metabolic_gene_edges", "core_gene_reaction",
    "reaction_membership", "meta_module_summary"
  )
  out <- lapply(names_to_merge, function(name) {
    .rc_bind_frames_fill(lapply(
      artifacts,
      function(artifact) artifact$grn_meta_modules[[name]]
    ))
  })
  names(out) <- names_to_merge

  core <- out$core_gene_reaction
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  core_ids <- unique(as.character(core$reaction_id))
  core_ids <- core_ids[!is.na(core_ids) & nzchar(core_ids)]

  biological <- out$reaction_membership
  biological_ids <- unique(as.character(biological$reaction_id))
  biological_ids <- biological_ids[
    !is.na(biological_ids) & nzchar(biological_ids)
  ]
  if (!length(core_ids) || !length(biological_ids)) {
    stop(
      "No merged biological meta-module reactions were produced.",
      call. = FALSE
    )
  }

  out$biological_reaction_membership <- biological
  out$merged_core_reactions <- data.frame(
    sample_id = "merged",
    module_id = "MERGED_META_MODULES",
    reaction_id = core_ids,
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  out$merged_reaction_membership <- data.frame(
    sample_id = "merged",
    module_id = "MERGED_META_MODULES",
    reaction_id = biological_ids,
    is_core = biological_ids %in% core_ids,
    inclusion_stage = ifelse(
      biological_ids %in% core_ids,
      "merged_meta_module_core",
      "merged_meta_module_biological_member"
    ),
    stringsAsFactors = FALSE
  )
  out$schema_version <- "regcompass_merged_meta_modules_v1"
  out$source_group_ids <- unique(vapply(
    artifacts,
    function(artifact) as.character(artifact$group_id),
    character(1)
  ))
  out$merge_source <- "deduplicated_biological_meta_module_reactions"
  out$is_gem <- FALSE
  out$fastcore_applied <- FALSE
  out
}
