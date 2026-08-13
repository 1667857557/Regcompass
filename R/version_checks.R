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

.rc_worker_runtime_profile <- function() {
  trimws(Sys.getenv("REGCOMPASS_WORKER_PROFILE", unset = ""))
}

.rc_skip_eager_seurat_stack <- function() {
  identical(.rc_worker_runtime_profile(), "layer2_numeric")
}

.onLoad <- function(libname, pkgname) {
  if (isTRUE(.rc_skip_eager_seurat_stack())) return(invisible(NULL))
  .rc_validate_seurat_stack_versions(.rc_seurat_stack_versions(), error = TRUE)
  invisible(NULL)
}
