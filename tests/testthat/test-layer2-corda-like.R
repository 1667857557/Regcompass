test_that("CORDA2 options preserve FASTCORE and match Python defaults", {
  defaults <- RegCompassR:::.rc_layer2_corda_options(list())
  expect_identical(defaults$model_completion, "fastcore")
  corda2 <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda2"
  ))
  expect_identical(corda2$model_completion, "corda")
  expect_identical(corda2$requested_model_completion, "corda2")
  expect_equal(corda2$penalty_factor, 100)
  expect_equal(corda2$cost_increase, 1.01)
  expect_equal(corda2$target_flux, 1)
  expect_equal(corda2$redundancies, 3L)
  expect_equal(corda2$support, 5L)
  expect_equal(corda2$upper_bound, 1e6)
  expect_identical(
    corda2$algorithm,
    "resendislab_python_CORDA2_corrected_redundant_path_assessment"
  )
  for (alias in c("corda", "corda_like")) {
    value <- RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = alias
    ))
    expect_identical(value$requested_model_completion, "corda2")
  }
  expect_error(
    RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = "corda2",
      corda2_cost_increase = 1
    )),
    "greater than 1"
  )
})

test_that("CORDA2 evidence is calculated within each cell type", {
  reactions <- paste0("R", 1:4)
  units <- paste0("U", 1:4)
  rna <- matrix(
    c(
      0.1, 0.2, 0.9, 1.0,
      0.2, 0.3, 0.8, 0.9,
      0.3, 0.4, 0.7, 0.8,
      0.4, 0.5, 0.6, 0.7
    ),
    nrow = 4,
    dimnames = list(reactions, units)
  )
  regulatory <- matrix(
    0.5, nrow = 4, ncol = 4,
    dimnames = list(reactions, units)
  )
  layer1 <- list(
    reaction_expression = rna,
    reaction_expression_rna_only = rna,
    reaction_regulatory_support_fraction = regulatory,
    unit_meta = data.frame(
      unit_id = units,
      cell_type = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    )
  )
  meta_modules <- list(workflow_params = list(celltype_col = "cell_type"))
  evidence <- RegCompassR:::.rc_layer2_corda_reaction_evidence(
    layer1, meta_modules, regulatory_weight = 0.2
  )
  expect_equal(nrow(evidence), 8L)
  expect_setequal(unique(evidence$cell_type), c("A", "B"))
  expect_true(all(evidence$evidence_score >= 0 &
                  evidence$evidence_score <= 1))
  expect_identical(
    unique(evidence$evidence_schema),
    "regcompass_corda_reaction_evidence_v2"
  )
})

test_that("RegCompass confidence mapping supplies five CORDA2 levels", {
  evidence <- data.frame(
    reaction_id = paste0("R", 1:7),
    evidence_score = c(0.9, 0.8, 0.7, 0.95, 0.4, 0.05, NA),
    stringsAsFactors = FALSE
  )
  classes <- RegCompassR:::.rc_corda_classify_reactions(
    parent_reactions = evidence$reaction_id,
    module_reactions = c("R1", "R2", "R3"),
    core_reactions = "R1",
    reaction_evidence = evidence,
    medium_confidence_threshold = 0.75,
    negative_confidence_threshold = 0.10,
    include_evidence_outside_modules = TRUE,
    max_medium_confidence_reactions = Inf
  )
  expect_identical(classes$hc, "R1")
  expect_setequal(classes$mc_module, c("R2", "R3"))
  expect_identical(classes$mc_evidence, "R4")
  expect_identical(classes$nc, "R6")
  expect_setequal(classes$ot, c("R5", "R7"))
  S <- Matrix::Diagonal(7)
  dimnames(S) <- list(paste0("M", 1:7), evidence$reaction_id)
  split <- RegCompassR:::.rc_corda_split_model(list(
    S = S,
    lb = stats::setNames(rep(0, 7), evidence$reaction_id),
    ub = stats::setNames(rep(10, 7), evidence$reaction_id)
  ))
  directional <- RegCompassR:::.rc_corda2_directional_confidence(
    split, classes
  )
  expect_equal(directional[["R1::forward"]], 3L)
  expect_equal(directional[["R2::forward"]], 2L)
  expect_equal(directional[["R4::forward"]], 1L)
  expect_equal(directional[["R6::forward"]], -1L)
  expect_equal(directional[["R7::forward"]], 0L)
})

test_that("direction splitting and CORDA2 bound normalization are explicit", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1, 1, -1), nrow = 2), sparse = TRUE,
    dimnames = list(c("A", "B"), c("IRR", "REV"))
  )
  gem <- list(
    S = S,
    lb = c(IRR = 0, REV = -5),
    ub = c(IRR = 10, REV = 7)
  )
  split <- RegCompassR:::.rc_corda_split_model(gem, tolerance = 1e-8)
  normalized <- RegCompassR:::.rc_corda2_normalize_split(split)
  expect_setequal(
    split$direction_table$variable_id,
    c("IRR::forward", "REV::forward", "REV::reverse")
  )
  expect_equal(normalized$ub[["IRR::forward"]], 1e6)
  expect_equal(normalized$ub[["REV::forward"]], 1e6)
  expect_equal(normalized$ub[["REV::reverse"]], 1e6)
})

test_that("Layer 2 exposes native CORDA2 without changing public formals", {
  implementation <- paste(
    deparse(body(rc_regcompass_step_layer2)), collapse = "\n"
  )
  expect_match(implementation, "corda2", fixed = TRUE)
  expect_match(
    implementation,
    ".rc_regcompass_step_layer2_corda_pool_base",
    fixed = TRUE
  )
  expect_identical(
    names(formals(rc_regcompass_step_layer2)),
    c(
      "layer1", "meta_modules", "gem", "medium_scenarios", "outdir",
      "model_mode", "layer2_args", "parallel", "BPPARAM", "progress"
    )
  )
})

test_that("CORDA2 cache is native and FASTCORE fallback is exact", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_build_celltype_medium_union_gem_cache)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    ".rc_build_celltype_medium_union_gem_cache_fastcore",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "cache_dir <- file.path(cache_dir, \"corda2\")",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "celltype_medium_corrected_python_corda2",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "attr(cache, \"completion_method\") <- \"corda2\"",
    fixed = TRUE
  )
  expect_false(grepl("file.rename", implementation, fixed = TRUE))
  expect_false(grepl("file.copy", implementation, fixed = TRUE))
})

test_that("CORDA2 runtime writes one checksum per shared structural file", {
  implementation <- paste(
    deparse(body(RegCompassR:::.rc_build_celltype_medium_union_gem_cache)),
    collapse = "\n"
  )
  expect_match(
    implementation,
    "checksum <- unname(tools::md5sum(file))",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "file_checksum = checksum",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "condition = \"all\"",
    fixed = TRUE
  )
})
