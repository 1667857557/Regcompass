.rc_with_preserved_seed <- function(seed, code) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || abs(seed - round(seed)) > sqrt(.Machine$double.eps)) {
    stop("`sample_balance_seed` must be one non-negative integer.", call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.rc_balance_condition_celltype_cells <- function(
    object, sample_col, condition_col, celltype_col,
    sample_balance = TRUE, sample_balance_seed = 12345L) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!is.logical(sample_balance) || length(sample_balance) != 1L ||
      is.na(sample_balance)) {
    stop("`sample_balance` must be TRUE or FALSE.", call. = FALSE)
  }
  required <- c(sample_col, condition_col, celltype_col)
  missing <- setdiff(required, colnames(object@meta.data))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  meta <- object@meta.data
  meta$.rc_cell_id <- rownames(meta)
  meta$.rc_sample <- trimws(as.character(meta[[sample_col]]))
  meta$.rc_condition <- trimws(as.character(meta[[condition_col]]))
  meta$.rc_celltype <- trimws(as.character(meta[[celltype_col]]))
  if (anyNA(meta[, c(".rc_sample", ".rc_condition", ".rc_celltype")]) ||
      any(!nzchar(meta$.rc_sample)) || any(!nzchar(meta$.rc_condition)) ||
      any(!nzchar(meta$.rc_celltype))) {
    stop("Sample, condition, and cell-type metadata must be complete.",
         call. = FALSE)
  }
  meta$.rc_balance_group <- paste(
    meta$.rc_condition,
    meta$.rc_celltype,
    sep = "\037"
  )
  group_ids <- sort(unique(meta$.rc_balance_group))

  balance_one_group <- function(group_id) {
    group_meta <- meta[meta$.rc_balance_group == group_id, , drop = FALSE]
    sample_ids <- sort(unique(group_meta$.rc_sample))
    counts <- table(factor(group_meta$.rc_sample, levels = sample_ids))
    target <- if (isTRUE(sample_balance)) min(as.integer(counts)) else NA_integer_
    retained <- lapply(sample_ids, function(sample_id) {
      cells <- sort(group_meta$.rc_cell_id[group_meta$.rc_sample == sample_id])
      if (!isTRUE(sample_balance) || length(cells) <= target) return(cells)
      cells[sample.int(length(cells), size = target, replace = FALSE)]
    })
    names(retained) <- sample_ids
    diagnostics <- do.call(rbind, lapply(sample_ids, function(sample_id) {
      n_input <- as.integer(counts[[sample_id]])
      n_retained <- length(retained[[sample_id]])
      data.frame(
        condition = group_meta$.rc_condition[[1L]],
        cell_type = group_meta$.rc_celltype[[1L]],
        biological_sample_id = sample_id,
        n_input_cells = n_input,
        target_cells_per_sample = if (isTRUE(sample_balance)) target else NA_integer_,
        n_retained_cells = n_retained,
        n_excluded_cells = n_input - n_retained,
        sample_balance = sample_balance,
        balance_strategy = if (isTRUE(sample_balance)) {
          "equal_cells_per_sample_at_minimum_sample_count"
        } else {
          "disabled_cell_count_weighted"
        },
        stringsAsFactors = FALSE
      )
    }))
    list(cells = unlist(retained, use.names = FALSE), diagnostics = diagnostics)
  }

  balanced <- .rc_with_preserved_seed(
    sample_balance_seed,
    lapply(group_ids, balance_one_group)
  )
  diagnostics <- do.call(rbind, lapply(balanced, `[[`, "diagnostics"))
  keep_cells <- unlist(lapply(balanced, `[[`, "cells"), use.names = FALSE)
  keep_cells <- keep_cells[!duplicated(keep_cells)]
  if (!length(keep_cells)) {
    stop("Sample balancing retained no cells.", call. = FALSE)
  }
  if (isTRUE(sample_balance) && any(diagnostics$n_retained_cells <= 0L)) {
    stop("Sample balancing produced an empty biological-sample stratum.",
         call. = FALSE)
  }
  balanced_object <- if (length(keep_cells) == ncol(object)) {
    object
  } else {
    subset(object, cells = keep_cells)
  }

  list(
    object = balanced_object,
    diagnostics = diagnostics,
    sample_balance = sample_balance,
    sample_balance_seed = as.integer(sample_balance_seed),
    sample_weighting = if (isTRUE(sample_balance)) {
      "equal_cells_per_sample_within_condition_celltype"
    } else {
      "cell_count_weighted"
    },
    n_input_cells = ncol(object),
    n_retained_cells = ncol(balanced_object),
    n_excluded_cells = ncol(object) - ncol(balanced_object)
  )
}

.rc_make_condition_pooled_metacells_unbalanced <-
  .rc_make_condition_pooled_metacells

.rc_make_condition_pooled_metacells_v170 <- function(
    object, outdir,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list(),
    strict_biological_defaults = TRUE) {
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
  }
  sample_balance <- metacell_args$sample_balance %||% TRUE
  sample_balance_seed <- metacell_args$sample_balance_seed %||% 12345L
  metacell_args$sample_balance <- NULL
  metacell_args$sample_balance_seed <- NULL

  balanced <- .rc_balance_condition_celltype_cells(
    object = object,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    sample_balance = sample_balance,
    sample_balance_seed = sample_balance_seed
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(
    balanced$diagnostics,
    file.path(outdir, "sample_balance_diagnostics.tsv.gz")
  )

  pooled <- .rc_make_condition_pooled_metacells_unbalanced(
    object = balanced$object,
    outdir = outdir,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args,
    strict_biological_defaults = strict_biological_defaults
  )
  pooled$sample_balance <- balanced$sample_balance
  pooled$sample_balance_seed <- balanced$sample_balance_seed
  pooled$sample_weighting <- balanced$sample_weighting
  pooled$sample_balance_diagnostics <- balanced$diagnostics
  pooled$sample_balance_summary <- balanced[c(
    "sample_balance", "sample_balance_seed", "sample_weighting",
    "n_input_cells", "n_retained_cells", "n_excluded_cells"
  )]
  pooled$metacell_meta$sample_weighting <- balanced$sample_weighting
  pooled$metacell_meta$sample_balance <- balanced$sample_balance
  pooled$input_design$sample_balance <- pooled$sample_balance_summary
  pooled
}

.rc_make_condition_pooled_metacells <-
  .rc_make_condition_pooled_metacells_v170
