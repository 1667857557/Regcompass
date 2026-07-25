.rc_seurat_object_export <- function(name) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) return(NULL)
  namespace <- asNamespace("SeuratObject")
  if (!exists(name, envir = namespace, inherits = FALSE)) return(NULL)
  get(name, envir = namespace, inherits = FALSE)
}

.rc_seurat_assay_names <- function(object) {
  assays <- tryCatch(
    SeuratObject::Assays(object),
    error = function(e) names(object@assays)
  )
  as.character(assays)
}

.rc_require_seurat_assay <- function(object, assay) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  assay <- trimws(as.character(assay[[1L]]))
  if (!nzchar(assay) || !assay %in% .rc_seurat_assay_names(object)) {
    stop("Seurat object is missing assay `", assay, "`.", call. = FALSE)
  }
  assay
}

.rc_assay_layer_names <- function(assay_object) {
  layers_fun <- .rc_seurat_object_export("Layers")
  if (!is.null(layers_fun)) {
    layers <- tryCatch(
      layers_fun(assay_object, search = NA),
      error = function(e) NULL
    )
    if (!is.null(layers)) return(as.character(layers))
  }
  intersect(
    c("counts", "data", "scale.data"),
    methods::slotNames(assay_object)
  )
}

.rc_matching_assay_layers <- function(assay_object, layer) {
  layers <- .rc_assay_layer_names(assay_object)
  layers[layers == layer | startsWith(layers, paste0(layer, "."))]
}

.rc_join_assay_layer <- function(object, assay, layer, required = TRUE) {
  assay <- .rc_require_seurat_assay(object, assay)
  assay_object <- object[[assay]]
  if (!inherits(assay_object, "Assay5")) {
    if (isTRUE(required) && !layer %in% .rc_assay_layer_names(assay_object)) {
      stop("Assay `", assay, "` has no `", layer, "` matrix.", call. = FALSE)
    }
    return(list(object = object, joined_layers = character()))
  }

  matches <- .rc_matching_assay_layers(assay_object, layer)
  if (!length(matches)) {
    if (!isTRUE(required)) {
      return(list(object = object, joined_layers = character()))
    }
    stop(
      "Assay `", assay, "` has no `", layer, "` layer.",
      call. = FALSE
    )
  }
  if (layer %in% matches && length(matches) > 1L) {
    stop(
      "Assay `", assay, "` contains both `", layer,
      "` and split `", layer, ".*` layers. Join or remove the ambiguous ",
      "layers before running RegCompass.",
      call. = FALSE
    )
  }
  if (identical(matches, layer)) {
    return(list(object = object, joined_layers = character()))
  }

  join_fun <- .rc_seurat_object_export("JoinLayers")
  if (is.null(join_fun)) {
    stop(
      "SeuratObject does not provide JoinLayers(), but assay `", assay,
      "` contains split `", layer, ".*` layers.",
      call. = FALSE
    )
  }
  object <- join_fun(
    object = object,
    assay = assay,
    layers = matches,
    new = layer
  )
  if (!layer %in% .rc_assay_layer_names(object[[assay]])) {
    stop(
      "Failed to join assay `", assay, "` layers into `", layer, "`.",
      call. = FALSE
    )
  }
  list(object = object, joined_layers = matches)
}

.rc_get_assay_matrix <- function(object, assay, layer) {
  assay <- .rc_require_seurat_assay(object, assay)
  assay_object <- object[[assay]]
  matches <- .rc_matching_assay_layers(assay_object, layer)
  if (!length(matches)) {
    stop("Assay `", assay, "` has no `", layer, "` matrix.", call. = FALSE)
  }
  if (length(matches) > 1L) {
    stop(
      "Assay `", assay, "` contains multiple `", layer,
      "` layers. RegCompass requires one joined layer; call JoinLayers() ",
      "or use the canonical RegCompass stage functions, which join a working ",
      "copy automatically.",
      call. = FALSE
    )
  }
  selected_layer <- matches[[1L]]

  layer_fun <- .rc_seurat_object_export("LayerData")
  value <- if (!is.null(layer_fun)) {
    layer_fun(object = object, assay = assay, layer = selected_layer)
  } else {
    SeuratObject::GetAssayData(
      object = object,
      assay = assay,
      slot = selected_layer
    )
  }
  if (is.null(dim(value)) || length(dim(value)) != 2L ||
      is.null(rownames(value)) || is.null(colnames(value))) {
    stop(
      "Assay `", assay, "` layer `", selected_layer,
      "` must be a named feature-by-cell matrix.",
      call. = FALSE
    )
  }
  value
}

.rc_set_assay_matrix <- function(object, assay, layer, new_data) {
  assay <- .rc_require_seurat_assay(object, assay)
  if (is.null(dim(new_data)) || length(dim(new_data)) != 2L ||
      is.null(rownames(new_data)) || is.null(colnames(new_data))) {
    stop("`new_data` must be a named feature-by-cell matrix.", call. = FALSE)
  }

  layer_setter <- .rc_seurat_object_export("LayerData<-")
  if (!is.null(layer_setter)) {
    return(layer_setter(
      object = object,
      assay = assay,
      layer = layer,
      value = new_data
    ))
  }
  SeuratObject::SetAssayData(
    object = object,
    assay = assay,
    slot = layer,
    new.data = new_data
  )
}

.rc_seurat_compatibility_summary <- function(
    object, assays, joined_layers = list()) {
  assays <- intersect(as.character(assays), .rc_seurat_assay_names(object))
  versions <- vapply(
    c("SeuratObject", "Seurat", "Signac"),
    function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(package))
    },
    character(1)
  )
  assay_classes <- stats::setNames(
    vapply(
      assays,
      function(assay) paste(class(object[[assay]]), collapse = "/"),
      character(1)
    ),
    assays
  )
  assay_storage <- stats::setNames(
    vapply(
      assays,
      function(assay) {
        if (inherits(object[[assay]], "Assay5")) "v5_layers" else "v3_slots"
      },
      character(1)
    ),
    assays
  )
  object_version <- tryCatch(
    as.character(SeuratObject::Version(object)),
    error = function(e) NA_character_
  )
  stack <- tryCatch(
    .rc_validate_seurat_stack_versions(versions, error = FALSE),
    error = function(e) list(profile = "unsupported", reason = conditionMessage(e))
  )
  list(
    default_input_profile = "Seurat_v4_v3_Assay",
    installed_stack_profile = stack$profile %||% "unsupported",
    package_versions = versions,
    object_version = object_version,
    assay_classes = assay_classes,
    assay_storage = assay_storage,
    joined_layers = joined_layers
  )
}

.rc_prepare_seurat_assays <- function(
    object, assays, required_layers = "counts", optional_layers = character()) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  assays <- unique(as.character(assays))
  missing <- setdiff(assays, .rc_seurat_assay_names(object))
  if (length(missing)) {
    stop("Seurat object is missing assays: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  previous <- tryCatch(
    object@misc$regcompass_seurat_compatibility,
    error = function(e) NULL
  )
  joined <- list()
  process_layer <- function(assay, layer, required) {
    result <- .rc_join_assay_layer(
      object, assay, layer, required = required
    )
    object <<- result$object
    if (length(result$joined_layers)) {
      joined[[paste(assay, layer, sep = ":")]] <<- result$joined_layers
    }
    if (isTRUE(required)) {
      invisible(.rc_get_assay_matrix(object, assay, layer))
    }
  }
  for (assay in assays) {
    for (layer in required_layers) process_layer(assay, layer, TRUE)
    for (layer in optional_layers) process_layer(assay, layer, FALSE)
  }
  previous_assays <- if (is.list(previous$assay_classes)) {
    names(previous$assay_classes)
  } else {
    names(previous$assay_classes %||% character())
  }
  all_assays <- unique(c(previous_assays, assays))
  previous_joined <- previous$joined_layers %||% list()
  all_joined <- c(previous_joined, joined)
  if (length(all_joined) && anyDuplicated(names(all_joined))) {
    all_joined <- all_joined[!duplicated(names(all_joined), fromLast = TRUE)]
  }
  object@misc$regcompass_seurat_compatibility <-
    .rc_seurat_compatibility_summary(object, all_assays, all_joined)
  object
}

.rc_get_assay_counts <- function(object, assay) {
  .rc_get_assay_matrix(object, assay, "counts")
}
