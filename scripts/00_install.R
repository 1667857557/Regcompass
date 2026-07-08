required <- c("Seurat", "SeuratObject", "testthat", "Matrix")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  install.packages(missing)
}

supercell_repo <- "1667857557/SuperCell-Seurat-V4"
supercell_ref <- "supercell-2.0"
supercell_url <- paste0("https://github.com/", supercell_repo, "/tree/", supercell_ref)

if (!requireNamespace("SuperCell", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github(supercell_repo, ref = supercell_ref, upgrade = "never")
}

message("Using SuperCell from ", supercell_url, " when SuperCell needs to be installed.")

if (utils::packageVersion("Seurat") >= "5.0.0") {
  warning("RegCompassR v0.1 is developed against Seurat v4-style objects and slots.")
}
