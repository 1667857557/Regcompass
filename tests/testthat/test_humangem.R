test_that("Human-GEM preparation returns symbol GPR tables", {
  repo <- file.path(tempdir(), paste0("Human-GEM-mini-", Sys.getpid()))
  model <- file.path(repo, "model")
  dir.create(model, recursive = TRUE)
  writeLines(c(
    "genes\tgeneSymbols",
    "ENSG00000156515\tHK1",
    "ENSG00000003436\tPFKM",
    "ENSG00000134333\tLDHA"
  ), file.path(model, "genes.tsv"))
  writeLines(c(
    "rxns\trxnKEGGID",
    "MAR00001\tR00001",
    "MAR00002\tR00002"
  ), file.path(model, "reactions.tsv"))
  writeLines(c(
    "reactions:",
    "- id: MAR00001",
    "  name: hexokinase",
    "  - gene_reaction_rule: ENSG00000156515 or ENSG00000003436",
    "- id: MAR00002",
    "  name: lactate dehydrogenase",
    "  - gene_reaction_rule: ENSG00000134333"
  ), file.path(model, "Human-GEM.yml"))

  out <- rc_prepare_humangem_gpr_table(repo)
  expect_true(all(c("gpr_table", "metabolic_genes", "reaction_rules", "genes", "reactions") %in% names(out)))
  expect_setequal(out$gpr_table$gene, c("HK1", "PFKM", "LDHA"))
  expect_setequal(out$metabolic_genes, c("HK1", "PFKM", "LDHA"))
  expect_true(all(c("reaction_id", "and_group_id", "gene") %in% colnames(out$gpr_table)))
})

test_that("Human-GEM preparation can keep Ensembl gene IDs", {
  repo <- file.path(tempdir(), paste0("Human-GEM-mini-ensembl-", Sys.getpid()))
  model <- file.path(repo, "model")
  dir.create(model, recursive = TRUE)
  writeLines(c("genes\tgeneSymbols", "ENSG1\tHK1"), file.path(model, "genes.tsv"))
  writeLines(c("rxns", "MAR00001"), file.path(model, "reactions.tsv"))
  writeLines(c("reactions:", "- id: MAR00001", "  - gene_reaction_rule: ENSG1"), file.path(model, "Human-GEM.yml"))

  out <- rc_prepare_humangem_gpr_table(repo, gene_format = "ensembl")
  expect_identical(out$gpr_table$gene, "ENSG1")
})

test_that("Human-GEM YAML parser supports omap/list-prefixed gene_reaction_rule", {
  yml <- tempfile(fileext = ".yml")
  writeLines(c(
    "- id: MAR00001",
    "  - gene_reaction_rule: ENSG000001 or ENSG000002"
  ), yml)
  out <- RegCompassR:::rc_read_humangem_yml_rules(yml)
  expect_equal(out$reaction_id, "MAR00001")
  expect_equal(out$gpr, "ENSG000001 or ENSG000002")
})
