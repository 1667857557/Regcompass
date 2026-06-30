required <- c("Seurat", "SeuratObject", "testthat", "Matrix")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  install.packages(missing)
}

if (utils::packageVersion("Seurat") >= "5.0.0") {
  warning("RegCompassR v0.1 is developed against Seurat v4-style objects and slots.")
}
