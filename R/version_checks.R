.rc_runtime_state <- new.env(parent = emptyenv())
.rc_runtime_state$namespace_loaded_at <- as.POSIXct(NA)

.rc_namespace_install_is_current <- function(
    loaded_at, installed_at, tolerance_seconds = 1) {
  loaded_at <- as.POSIXct(loaded_at)
  installed_at <- as.POSIXct(installed_at)
  tolerance_seconds <- suppressWarnings(as.numeric(tolerance_seconds))
  if (length(loaded_at) != 1L || is.na(loaded_at) ||
      length(installed_at) != 1L || is.na(installed_at) ||
      length(tolerance_seconds) != 1L || !is.finite(tolerance_seconds) ||
      tolerance_seconds < 0) {
    return(TRUE)
  }
  !isTRUE(installed_at > loaded_at + tolerance_seconds)
}

.rc_assert_current_namespace_install <- function() {
  installed_rdb <- system.file(
    "R", "RegCompassR.rdb", package = "RegCompassR"
  )
  if (!nzchar(installed_rdb) || !file.exists(installed_rdb)) {
    return(invisible(TRUE))
  }
  installed_at <- file.info(installed_rdb)$mtime
  loaded_at <- .rc_runtime_state$namespace_loaded_at
  if (!.rc_namespace_install_is_current(loaded_at, installed_at)) {
    stop(
      "RegCompassR was reinstalled after this R session loaded its namespace. ",
      "Restart R before running a RegCompass stage; otherwise the session ",
      "continues executing stale in-memory code.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_seurat_stack_versions <- function() {
  vapply(
    c("SeuratObject", "Seurat", "Signac"),
    function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(package))
    },
    character(1)
  )
}

.rc_validate_seurat_stack_versions <- function(observed, error = TRUE) {
  packages <- c("SeuratObject", "Seurat", "Signac")
  observed <- observed[packages]
  missing <- packages[is.na(observed) | !nzchar(observed)]
  fail <- function(message) {
    if (isTRUE(error)) stop(message, call. = FALSE)
    list(profile = "unsupported", supported = FALSE, reason = message)
  }
  if (length(missing)) {
    return(fail(paste0(
      "RegCompassR requires installed packages: ",
      paste(missing, collapse = ", "), "."
    )))
  }

  versions <- lapply(observed, numeric_version)
  names(versions) <- packages
  major <- function(value) suppressWarnings(
    as.integer(sub("[.].*$", "", as.character(value)))
  )
  seurat_object_major <- major(observed[["SeuratObject"]])
  seurat_major <- major(observed[["Seurat"]])
  signac_major <- major(observed[["Signac"]])

  if (!identical(seurat_object_major, seurat_major)) {
    return(fail(paste0(
      "SeuratObject and Seurat must use the same major version. Observed: ",
      paste(names(observed), observed, sep = "=", collapse = ", "), "."
    )))
  }
  if (!identical(signac_major, 1L)) {
    return(fail(paste0(
      "RegCompassR supports Signac 1.x ChromatinAssay objects; Signac 2.x ",
      "ChromatinAssay5 is not yet supported. Observed Signac=", observed[["Signac"]],
      "."
    )))
  }

  if (identical(seurat_major, 4L)) {
    valid <- versions$SeuratObject >= numeric_version("4.1.4") &&
      versions$SeuratObject < numeric_version("5.0.0") &&
      versions$Seurat >= numeric_version("4.4.0") &&
      versions$Seurat < numeric_version("5.0.0") &&
      versions$Signac >= numeric_version("1.11.0") &&
      versions$Signac < numeric_version("1.14.0")
    if (!valid) {
      return(fail(paste0(
        "Unsupported Seurat v4 stack. RegCompassR requires SeuratObject ",
        ">=4.1.4,<5; Seurat >=4.4.0,<5; and Signac >=1.11.0,<1.14.0. ",
        "Observed: ",
        paste(names(observed), observed, sep = "=", collapse = ", "), "."
      )))
    }
    return(list(
      profile = "seurat_v4_default",
      supported = TRUE,
      default = TRUE,
      observed = observed,
      assay_policy = "v3_Assay_and_ChromatinAssay"
    ))
  }

  if (identical(seurat_major, 5L)) {
    valid <- versions$SeuratObject >= numeric_version("5.0.0") &&
      versions$SeuratObject < numeric_version("6.0.0") &&
      versions$Seurat >= numeric_version("5.0.0") &&
      versions$Seurat < numeric_version("6.0.0") &&
      versions$Signac >= numeric_version("1.12.0") &&
      versions$Signac < numeric_version("2.0.0")
    if (!valid) {
      return(fail(paste0(
        "Unsupported Seurat v5 stack. RegCompassR requires SeuratObject ",
        ">=5,<6; Seurat >=5,<6; and Signac >=1.12.0,<2. Observed: ",
        paste(names(observed), observed, sep = "=", collapse = ", "), "."
      )))
    }
    return(list(
      profile = "seurat_v5_compatible",
      supported = TRUE,
      default = FALSE,
      observed = observed,
      assay_policy = "v3_Assay_or_joinable_Assay5_plus_ChromatinAssay"
    ))
  }

  fail(paste0(
    "RegCompassR supports Seurat major versions 4 and 5 only. Observed: ",
    paste(names(observed), observed, sep = "=", collapse = ", "), "."
  ))
}

.onLoad <- function(libname, pkgname) {
  .rc_runtime_state$namespace_loaded_at <- Sys.time()
  .rc_validate_seurat_stack_versions(.rc_seurat_stack_versions(), error = TRUE)
}
