test_that("metacell import preserves single-letter metadata labels", {
  directory <- tempfile("regcompass-metacell-")
  dir.create(directory, recursive = TRUE)

  metadata <- data.frame(
    metacell_id = c("MC1", "MC2"),
    sample_id = c("S1", "S1"),
    condition = c("T", "F"),
    cell_type = c("T", "T"),
    n_cells = c(10L, 12L),
    stringsAsFactors = FALSE
  )
  write_gz <- function(x, path) {
    connection <- gzfile(path, open = "wt")
    on.exit(close(connection), add = TRUE)
    utils::write.table(
      x,
      file = connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  write_gz(metadata, file.path(directory, "metacell_metadata.tsv.gz"))

  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3),
      nrow = 2,
      dimnames = list(c("G1", "G2"), c("MC1", "MC2"))
    ),
    sparse = TRUE
  )
  saveRDS(counts, file.path(directory, "rna_counts.rds"))

  imported <- rc_import_supercell2_metacells(directory)

  expect_identical(imported$metacell_meta$condition, c("T", "F"))
  expect_identical(imported$metacell_meta$cell_type, c("T", "T"))
  expect_false(is.logical(imported$metacell_meta$condition))
  expect_equal(imported$metacell_meta$n_cells, c(10L, 12L))
})

test_that("database source annotations omit empty prefixes", {
  S <- diag(4)
  dimnames(S) <- list(paste0("M", 1:4), paste0("R", 1:4))
  reaction_meta <- data.frame(
    reaction_id = paste0("R", 1:4),
    subsystem = c("A", "A", "B", "C"),
    metabolic_module = c("A", "A", "B", "C"),
    kegg_reaction_id = c("K1", NA, "K1", NA),
    reactome_reaction_id = c(NA, "X1", NA, "X1"),
    rhea_master_id = NA_character_,
    stringsAsFactors = FALSE
  )
  gem <- rc_make_gem(
    S,
    lb = rep(0, 4),
    ub = rep(1000, 4),
    reaction_meta = reaction_meta
  )
  core <- data.frame(
    sample_id = "S1",
    module_id = "S1::GRN0001",
    gene = "G1",
    reaction_id = "R1",
    stringsAsFactors = FALSE
  )

  expanded <- rc_expand_meta_module_reactions(gem, core)
  source <- stats::setNames(
    expanded$reaction_membership$source_annotation,
    expanded$reaction_membership$reaction_id
  )

  expect_identical(source[["R3"]], "KEGG:K1")
  expect_identical(source[["R4"]], "REACTOME:X1")
  expect_false(any(grepl(
    "(^|;)(KEGG|REACTOME):(;|$)",
    stats::na.omit(expanded$reaction_membership$source_annotation)
  )))
})
