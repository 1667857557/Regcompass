rc_load_or_merge_metacell_objects <- function(
    objects, fragment_manifest = NULL, metacell_meta = NULL,
    rna_assay = "RNA", atac_assay = "ATAC",
    require_complete_fragments = FALSE) {
  if (!is.list(objects) && !is.character(objects)) {
    stop("`objects` must contain Seurat objects or RDS paths.", call. = FALSE)
  }
  values <- lapply(objects, function(value) {
    if (inherits(value, "Seurat")) return(value)
    if (is.character(value) && length(value) == 1L && file.exists(value)) {
      object <- readRDS(value)
      if (!inherits(object, "Seurat")) {
        stop("Metacell RDS does not contain a Seurat object: ", value,
             call. = FALSE)
      }
      return(object)
    }
    stop("Invalid metacell object entry.", call. = FALSE)
  })
  if (!length(values)) stop("No metacell objects were supplied.", call. = FALSE)

  all_ids <- unlist(lapply(values, colnames), use.names = FALSE)
  duplicated_ids <- unique(all_ids[duplicated(all_ids)])
  if (length(duplicated_ids)) {
    stop(
      "Metacell IDs are not globally unique before Seurat merge: ",
      paste(utils::head(duplicated_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  removed_extra <- character()
  expected_ids <- NULL
  id_col <- NULL
  if (!is.null(metacell_meta)) {
    if (!is.data.frame(metacell_meta)) {
      stop("`metacell_meta` must be a data frame.", call. = FALSE)
    }
    id_col <- if ("metacell_id" %in% colnames(metacell_meta)) {
      "metacell_id"
    } else if ("pool_id" %in% colnames(metacell_meta)) {
      "pool_id"
    } else {
      stop("Metacell metadata require metacell_id or pool_id.", call. = FALSE)
    }
    expected_ids <- as.character(metacell_meta[[id_col]])
    if (anyNA(expected_ids) || any(!nzchar(trimws(expected_ids))) ||
        anyDuplicated(expected_ids)) {
      stop("Expected metacell IDs must be unique and non-empty.",
           call. = FALSE)
    }
    missing <- setdiff(expected_ids, all_ids)
    if (length(missing)) {
      stop(
        "Expected metacell IDs are absent from supplied objects: ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    removed_extra <- setdiff(all_ids, expected_ids)
    values <- lapply(values, function(value) {
      keep <- intersect(colnames(value), expected_ids)
      if (!length(keep)) return(NULL)
      subset(value, cells = keep)
    })
    values <- values[!vapply(values, is.null, logical(1))]
  }

  object <- if (length(values) == 1L) {
    values[[1L]]
  } else {
    Reduce(function(x, y) merge(x, y, merge.data = TRUE), values)
  }

  if (!is.null(expected_ids)) {
    object <- subset(object, cells = expected_ids)
    index <- match(expected_ids, as.character(metacell_meta[[id_col]]))
    object@meta.data <- metacell_meta[index, , drop = FALSE]
    rownames(object@meta.data) <- expected_ids
  }
  attr(object, "removed_extra_metacell_ids") <- removed_extra

  if (isTRUE(require_complete_fragments) && is.null(fragment_manifest)) {
    stop("Complete fragment registration was required but no manifest was supplied.",
         call. = FALSE)
  }
  if (!is.null(fragment_manifest)) {
    registration <- .rc_fragment_registration_from_manifest(
      fragment_manifest, object_cells = colnames(object)
    )
    object@misc$regcompass_fragment_registration <- registration
  }
  object
}
