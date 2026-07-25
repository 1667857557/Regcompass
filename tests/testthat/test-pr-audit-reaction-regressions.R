audit_reaction_fixture <- function() {
  S <- matrix(
    c(
      -1, 0,
      -2, 0,
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

  conditions <- rep(c("control", "MS177", "JQ1"), each = 6L)
  units <- paste0("audit_u", seq_along(conditions))
  rna <- matrix(
    0.5,
    nrow = 2,
    ncol = length(units),
    dimnames = list(c("GENE1", "GENE2"), units)
  )
  modifier <- matrix(0, nrow = 2, ncol = length(units), dimnames = dimnames(rna))
  modifier["GENE1", conditions == "MS177"] <- 0.5
  modifier["GENE1", conditions == "JQ1"] <- 0.2
  multiome <- rna
  multiome["GENE1", conditions == "MS177"] <- 0.65
  multiome["GENE1", conditions == "JQ1"] <- 0.56
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
    parsed_gpr = list(R_GENE = list(c("GENE1", "GENE2"))),
    capacity_params = list(
      promiscuity_mode = "none",
      tau = 0.20,
      and_method = "boltzmann",
      or_method = "sum"
    )
  )
  row_id <- "reaction=R_GENE::direction=forward::medium=base"
  penalty_values <- c(seq(12, 17), seq(1, 6), seq(7, 12))
  penalty <- matrix(
    penalty_values,
    nrow = 1,
    dimnames = list(row_id, units)
  )
  microcompass <- list(
    penalty = penalty,
    vmax = matrix(
      100, nrow = 1, ncol = length(units), dimnames = dimnames(penalty)
    ),
    feasible = matrix(
      TRUE, nrow = 1, ncol = length(units), dimnames = dimnames(penalty)
    ),
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
  result <- rc_attach_reaction_annotations(
    result,
    gem,
    condition_col = "condition",
    celltype_col = "cell_type"
  )
  list(gem = gem, layer1 = layer1, result = result, row_id = row_id)
}

test_that("gene-level fallback cannot claim active RNA+ATAC reaction evidence", {
  units <- c("u1", "u2")
  layer1 <- list(
    gene_support_rna = matrix(
      c(1, 1), nrow = 1, dimnames = list("Slc22a17", units)
    ),
    gene_regulatory_modifier = matrix(
      c(0.5, 0.5), nrow = 1, dimnames = list("Slc22a17", units)
    ),
    gene_support_multiome = matrix(
      c(1.5, 1.5), nrow = 1, dimnames = list("Slc22a17", units)
    ),
    unit_meta = data.frame(
      pool_id = units,
      condition = "A",
      cell_type = "C",
      stringsAsFactors = FALSE
    )
  )
  catalog <- data.frame(
    reaction_id = "R1",
    genes = "Slc22a17",
    stringsAsFactors = FALSE
  )
  evidence <- .rc_ra_group_evidence(
    catalog,
    layer1,
    condition_col = "condition",
    celltype_col = "cell_type"
  )
  expect_identical(evidence$evidence_class, "RNA-only")
  expect_identical(
    evidence$evidence_resolution,
    "reaction_capacity_unavailable"
  )
  expect_false(evidence$has_active_multiome_contribution)
  expect_identical(evidence$atac_modifier_genes, "Slc22a17")
  expect_identical(evidence$multiome_contributing_genes, "Slc22a17")

  omnibus <- .rc_ra_omnibus_evidence(
    data.frame(reaction_id = "missing", cell_type = "C"),
    evidence,
    conditions = c("A", "B")
  )
  expect_identical(omnibus$evidence_class, "unknown/unavailable")
})

test_that("reaction catalogs use normalized bounds inferred roles and source gene case", {
  S <- matrix(
    c(
      -1, 0,
      0, -1,
      0, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A_e", "B_c", "C_c"), c("EX_A", "R_MOUSE"))
  )
  gem <- list(
    S = S,
    lb = c(-10, 0),
    ub = c(1000, 1000),
    reaction_meta = data.frame(
      reaction_id = c("EX_A", "R_MOUSE"),
      name = c("A exchange", "Mouse conversion"),
      stringsAsFactors = FALSE
    ),
    metabolite_meta = data.frame(
      metabolite_id = c("A_e", "B_c", "C_c"),
      name = c("A", "B", "C"),
      compartment = c("e", "c", "c"),
      stringsAsFactors = FALSE
    ),
    gpr_table = data.frame(
      reaction_id = "R_MOUSE",
      and_group_id = 1L,
      gene = "Slc22a17",
      stringsAsFactors = FALSE
    )
  )
  catalog <- .rc_ra_reaction_catalog(gem)
  expect_identical(
    catalog$reaction_role[catalog$reaction_id == "EX_A"],
    "exchange"
  )
  expect_equal(
    catalog$lower_bound[catalog$reaction_id == "R_MOUSE"],
    0
  )
  expect_identical(
    catalog$genes[catalog$reaction_id == "R_MOUSE"],
    "Slc22a17"
  )
  selected <- rc_select_gene_reactions(
    list(reaction_catalog = catalog, reaction_evidence = data.frame()),
    genes = "SLC22A17"
  )
  expect_identical(selected$reactions$matched_genes, "Slc22a17")
})

test_that("condition plots preserve reaction annotation context", {
  skip_if_not_installed("ggplot2")
  fixture <- audit_reaction_fixture()
  plot <- rc_plot_condition_reaction(
    fixture$result,
    reaction_id = "R_GENE",
    cell_type = "stem-cell_like",
    target_direction = "forward",
    medium_scenario = "base",
    condition_col = "condition",
    celltype_col = "cell_type",
    conditions = c("control", "MS177"),
    min_units = 5L
  )
  statistics <- attr(plot, "condition_statistics")
  annotation <- attr(plot, "reaction_annotation")
  expect_true("reaction_name" %in% colnames(statistics$pairwise))
  expect_identical(annotation$reaction_name, "Gene-associated conversion")
  expect_match(annotation$evidence_comparison, "RNA-only")
  expect_match(annotation$evidence_comparison, "RNA\\+ATAC")
})

test_that("gene plot filters use requested conditions and min_units", {
  skip_if_not_installed("ggplot2")
  fixture <- audit_reaction_fixture()
  plots <- rc_plot_condition_gene_reactions(
    fixture$result,
    genes = "GENE1",
    cell_type = "stem-cell_like",
    condition_col = "condition",
    celltype_col = "cell_type",
    conditions = c("control", "MS177"),
    min_units = 5L,
    target_directions = "forward",
    medium_scenario = "base",
    p_adj_max = 0.05,
    min_abs_rank_biserial = 0.3
  )
  expect_setequal(
    unique(plots$gene_selection$evidence$condition),
    c("control", "MS177")
  )
  expect_identical(plots$statistics$params$min_units, 5L)
  expect_error(
    rc_plot_condition_gene_reactions(
      fixture$result,
      genes = "GENE1",
      cell_type = "stem-cell_like",
      condition_col = "condition",
      celltype_col = "cell_type",
      conditions = c("control", "MS177"),
      min_units = 7L,
      target_directions = "forward",
      medium_scenario = "base",
      p_adj_max = 0.05,
      min_abs_rank_biserial = 0.3
    ),
    "No gene-associated targets passed"
  )
})
