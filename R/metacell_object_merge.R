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
  object <- if (length(values) == 1L) {
    values[[1L]]
  } else {
    Reduce(function(x, y) merge(x, y), values)
  }
  if (!is.null(metacell_meta) && is.data.frame(metacell_meta)) {
    id_col <- if ("metacell_id" %in% colnames(metacell_meta)) {
      "metacell_id"
    } else if ("pool_id" %in% colnames(metacell_meta)) {
      "pool_id"
    } else {
      NULL
    }
    if (!is.null(id_col)) {
      index <- match(colnames(object), as.character(metacell_meta[[id_col]]))
      if (anyNA(index)) {
        stop("Metacell metadata do not cover the merged object.",
             call. = FALSE)
      }
      object@meta.data <- metacell_meta[index, , drop = FALSE]
      rownames(object@meta.data) <- colnames(object)
    }
  }
  object
}
