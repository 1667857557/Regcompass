test_that("symmetric condition coding has a zero-sum deviation basis", {
  contrast <- RegCompassR:::.rc_multitask_contrast(
    c("treated", "control", "treated", "other")
  )

  expect_equal(rownames(contrast), c("control", "other", "treated"))
  expect_equal(rowSums(contrast), rep(0, 3), tolerance = 1e-12)
  expect_equal(colSums(contrast), rep(0, 3), tolerance = 1e-12)
  expect_equal(unname(contrast), diag(3) - matrix(1 / 3, 3, 3))
})

test_that("decoded condition deviations sum to zero and preserve global mean", {
  condition <- c("A", "B", "C")
  contrast <- RegCompassR:::.rc_multitask_contrast(condition)
  beta <- c(0.5, -0.25)
  gamma <- matrix(
    c(1, 2, -1, 0.5, 0.2, -0.4),
    nrow = 2,
    ncol = 3
  )
  decoded <- RegCompassR:::.rc_decode_multitask_coefficients(
    c(beta, as.numeric(gamma)),
    n_edges = 2,
    contrast = contrast
  )

  expect_equal(colSums(decoded$delta), c(0, 0), tolerance = 1e-12)
  expect_equal(colMeans(decoded$theta), beta, tolerance = 1e-12)
  expect_equal(decoded$theta, sweep(decoded$delta, 2, beta, "+"))
})

test_that("condition balancing gives every condition equal total loss weight", {
  condition <- c(rep("A", 10), rep("B", 25), rep("C", 5))
  weight <- RegCompassR:::.rc_condition_balanced_weights(condition)
  total <- tapply(weight, condition, sum)

  expect_equal(unname(total), rep(total[[1]], 3), tolerance = 1e-12)
  expect_equal(mean(weight), 1, tolerance = 1e-12)
})

test_that("multitask validation requires a ridge component", {
  expect_error(
    RegCompassR:::.rc_validate_multitask_grn_args(list(alpha = 1)),
    "non-zero ridge component"
  )
  expect_equal(
    RegCompassR:::.rc_validate_multitask_grn_args(list(alpha = 0.5))$alpha,
    0.5
  )
})

test_that("zero regulatory evidence exactly returns RNA-only support", {
  rna <- matrix(
    c(0, 0.2, 0.5, 1),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("u1", "u2"))
  )
  modifier <- matrix(0, nrow = 2, ncol = 2, dimnames = dimnames(rna))
  integrated <- RegCompassR:::.rc_integrate_regulatory_support(
    rna, modifier, alpha = 2
  )

  expect_equal(integrated, rna, tolerance = 1e-12)
})

test_that("different condition gene sets produce different complete-GPR cores", {
  genes <- data.frame(
    sample_id = c("A", "A", "A", "B"),
    module_id = c("MA", "MA", "MA", "MB"),
    gene = c("G1", "G2", "G3", "G1"),
    stringsAsFactors = FALSE
  )
  gpr <- data.frame(
    reaction_id = c("R_complex", "R_complex", "R_iso", "R_iso"),
    and_group_id = c("1", "1", "1", "2"),
    gene = c("G1", "G2", "G3", "G1"),
    stringsAsFactors = FALSE
  )
  mapped <- RegCompassR:::rc_map_meta_module_core_reactions(genes, gpr)

  a_core <- unique(mapped$reaction_id[
    mapped$sample_id == "A" & mapped$is_core %in% TRUE
  ])
  b_core <- unique(mapped$reaction_id[
    mapped$sample_id == "B" & mapped$is_core %in% TRUE
  ])
  expect_setequal(a_core, c("R_complex", "R_iso"))
  expect_setequal(b_core, "R_iso")
  expect_false("R_complex" %in% b_core)
})

test_that("merged catalogue preserves multitask provenance and one reaction union", {
  condition_modules <- list(
    grn_mode = "multitask_shared_backbone",
    celltype_fit_status = data.frame(cell_type = "T", status = "ok"),
    sample_status = data.frame(group_id = c("A_T", "B_T")),
    tf_peak_gene_candidates = data.frame(
      edge_universe_id = "u1", edge_id = "e1"
    ),
    tf_peak_gene_global = data.frame(edge_id = "e1"),
    tf_peak_gene_condition_all = data.frame(edge_id = "e1"),
    tf_peak_gene_all = data.frame(edge_id = "e1"),
    tf_peak_gene_significant = data.frame(edge_id = "e1"),
    condition_target_genes = data.frame(target = "G1"),
    target_model_diagnostics = data.frame(target = "G1"),
    stability_diagnostics = data.frame(edge_id = "e1"),
    supported_metabolic_genes = data.frame(gene = "G1"),
    core_gene_reaction = data.frame(
      reaction_id = c("R1", "R2"), is_core = TRUE
    ),
    reaction_membership = data.frame(
      reaction_id = c("R1", "R2", "R3")
    ),
    meta_module_summary = data.frame(module_id = "M")
  )
  merged <- RegCompassR:::.rc_merge_meta_module_catalogue(condition_modules)

  expect_identical(
    merged$schema_version,
    "regcompass_merged_multitask_meta_modules_v1"
  )
  expect_setequal(
    merged$merged_core_reactions$reaction_id,
    c("R1", "R2")
  )
  expect_setequal(
    merged$merged_reaction_membership$reaction_id,
    c("R1", "R2", "R3")
  )
  expect_identical(merged$source_edge_universe_ids, "u1")
})
