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
