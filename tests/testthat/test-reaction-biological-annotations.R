make_reaction_annotation_fixture <- function() {
  S <- matrix(
    c(
      -1,  0,
      -2,  0,
       1, -1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("A_c", "B_c", "C_e"),
      c("R_GENE", "EX_C")
    )
  )
  gem <- rc_make_gem(
    S,
    lb = c(R_GENE = -1000, EX_C = 0),
    ub = c(R_GENE = 1000, EX_C = 1000),
    reaction_meta = data.frame(
      reaction_id = c("R_GENE", "EX_C"),
      name = c("Gene-associated conversion", "C exchange"),
      subsystem = c("Test metabolism", "Exchange/demand reaction"),
      role = c("internal", "exchange"),
      stringsAsFactors = FALSE
    ),
    metabolite_meta = data.frame(
      metabolite_id = c("A_c", "B_c", "C_e"),
      name = c("metabolite A", "metabolite B", "metabolite C"),
      compartment = c("c", "c", "e"),
      stringsAsFactors = FALSE
    )
  )
  gem$gpr_table <- data.frame(
    reaction_id = c("R_GENE", "R_GENE"),
    and_group_id = c(1L, 1L),
    gene = c("GENE1", "GENE2"),
    stringsAsFactors = FALSE
  )

  conditions <- rep(c("control", "MS177"), each = 6L)
  units <- paste0("u", seq_along(conditions))
  rna <- matrix(
    0.5,
    nrow = 2,
    ncol = length(units),
    dimnames = list(c("gene1", "gene2"), units)
  )
  modifier <- matrix(
    0,
    nrow = 2,
    ncol = length(units),
    dimnames = dimnames(rna)
  )
  modifier["gene1", conditions == "MS177"] <- 0.5
  multiome <- rna
  multiome["gene1", conditions == "MS177"] <- 0.65
  layer1 <- list(
    gene_support_rna = rna,
    gene_regulatory_modifier = modifier,
    gene_support_multiome = multiome,
    unit_meta = data.frame(
      pool_id = units,
      condition = conditions,
      cell_type = "stem-cell_like",
      stringsAsFactors = FALSE
    ),
    parsed_gpr = list(
      R_GENE = list(c("GENE1", "GENE2"))
    )
  )

  row_id <- "reaction=R_GENE::direction=forward::medium=base"
  penalty <- matrix(
    c(seq(12, 17), seq(1, 6)),
    nrow = 1,
    dimnames = list(row_id, units)
  )
  microcompass <- list(
    penalty = penalty,
    vmax = matrix(100, nrow = 1, ncol = length(units), dimnames = dimnames(penalty)),
    feasible = matrix(TRUE, nrow = 1, ncol = length(units), dimnames = dimnames(penalty)),
    unit_meta = layer1$unit_meta,
    params = list(omega = 0.95, unit = "metacell")
  )
  result <- list(
    layer1 = layer1,
    microcompass = microcompass,
    reaction_ranking = data.frame(),
    condition_summary = data.frame(),
    condition_contrast = data.frame(),
    params = list()
  )
  list(gem = gem, layer1 = layer1, result = result, row_id = row_id)
}

test_that("reaction annotations contain names formulas GPRs and evidence classes", {
  fixture <- make_reaction_annotation_fixture()
  annotation <- rc_build_reaction_annotations(
    fixture$gem,
    fixture$layer1,
    condition_col = "condition",
    celltype_col = "cell_type"
  )

  gene_reaction <- annotation$reactions[
    annotation$reactions$reaction_id == "R_GENE", , drop = FALSE
  ]
  expect_equal(gene_reaction$reaction_name, "Gene-associated conversion")
  expect_equal(
    gene_reaction$forward_formula,
    "metabolite A [c] + 2 metabolite B [c] -> metabolite C [e]"
  )
  expect_equal(
    gene_reaction$reverse_formula,
    "metabolite C [e] -> metabolite A [c] + 2 metabolite B [c]"
  )
  expect_equal(gene_reaction$genes, "GENE1;GENE2")
  expect_equal(gene_reaction$gpr_rule, "(GENE1 and GENE2)")

  control <- subset(
    annotation$evidence,
    reaction_id == "R_GENE" & condition == "control"
  )
  ms177 <- subset(
    annotation$evidence,
    reaction_id == "R_GENE" & condition == "MS177"
  )
  structural <- subset(
    annotation$evidence,
    reaction_id == "EX_C" & condition == "control"
  )
  expect_equal(control$evidence_class, "RNA-only")
  expect_equal(ms177$evidence_class, "RNA+ATAC")
  expect_equal(ms177$multiome_contributing_genes, "GENE1")
  expect_true(ms177$has_active_multiome_contribution)
  expect_equal(structural$evidence_class, "structural/no-GPR")
})

test_that("condition statistics are enriched with reaction biology", {
  fixture <- make_reaction_annotation_fixture()
  result <- rc_attach_reaction_annotations(
    fixture$result,
    fixture$gem,
    condition_col = "condition",
    celltype_col = "cell_type"
  )
  statistics <- rc_test_condition_reactions(
    result,
    condition_col = "condition",
    celltype_col = "cell_type",
    conditions = c("control", "MS177"),
    cell_types = "stem-cell_like",
    min_units = 5L
  )
  expect_true(all(c(
    "reaction_name", "tested_formula", "genes", "gpr_rule",
    "evidence_class_a", "evidence_class_b", "evidence_comparison"
  ) %in% colnames(statistics$pairwise)))
  expect_equal(statistics$pairwise$reaction_name, "Gene-associated conversion")
  expect_equal(statistics$pairwise$evidence_class_a, "RNA-only")
  expect_equal(statistics$pairwise$evidence_class_b, "RNA+ATAC")
  expect_equal(
    statistics$pairwise$tested_formula,
    "metabolite A [c] + 2 metabolite B [c] -> metabolite C [e]"
  )
})

test_that("metabolic genes select reactions and significant plot collections", {
  fixture <- make_reaction_annotation_fixture()
  result <- rc_attach_reaction_annotations(
    fixture$result,
    fixture$gem,
    condition_col = "condition",
    celltype_col = "cell_type"
  )
  selected <- rc_select_gene_reactions(
    result,
    genes = "GENE1",
    cell_types = "stem-cell_like"
  )
  expect_equal(selected$reaction_ids, "R_GENE")
  expect_equal(selected$reactions$matched_genes, "GENE1")

  skip_if_not_installed("ggplot2")
  plots <- rc_plot_condition_gene_reactions(
    result,
    genes = "GENE1",
    cell_type = "stem-cell_like",
    condition_col = "condition",
    celltype_col = "cell_type",
    conditions = c("control", "MS177"),
    target_directions = "forward",
    medium_scenario = "base",
    p_adj_max = 0.05,
    min_abs_rank_biserial = 0.3,
    max_reactions = 3L
  )
  expect_s3_class(plots, "regcompass_gene_reaction_plots")
  expect_equal(length(plots$plots), 1L)
  expect_s3_class(plots$plots[[1L]], "ggplot")
  expect_equal(plots$selected_targets$reaction_id, "R_GENE")
})
