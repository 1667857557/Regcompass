.onLoad <- function(libname, pkgname) {
  expected <- c(
    SeuratObject = "4.1.4",
    Seurat = "4.4.0",
    Signac = "1.11.0"
  )
  observed <- vapply(
    names(expected),
    function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(package))
    },
    character(1)
  )
  mismatch <- is.na(observed) | observed != expected
  if (any(mismatch)) {
    stop(
      paste0(
        "RegCompassR requires the validated Seurat stack: ",
        paste(names(expected), expected, sep = "=", collapse = ", "),
        ". Observed: ",
        paste(names(observed), observed, sep = "=", collapse = ", "),
        ". Install the pinned versions before loading RegCompassR."
      ),
      call. = FALSE
    )
  }
}
