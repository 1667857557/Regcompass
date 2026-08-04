test_that("CORDA-like options preserve FASTCORE by default", {
  options <- RegCompassR:::.rc_layer2_corda_options(list())
  expect_identical(options$model_completion, "fastcore")
  expect_equal(options$medium_confidence_threshold, 0.75)
  expect_equal(options$negative_confidence_threshold, 0.10)
  expect_equal(options$other_penalty, 1)
  expect_equal(options$negative_penalty, 10)
  expect_error(
    RegCompassR:::.rc_layer2_corda_options(list(
      corda_other_penalty = 2,
      corda_negative_penalty = 1
    )),
    "greater than or equal"
  )
})

test_that("CORDA-like evidence is calculated within each cell type", {
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
  multiome <- rna
  regulatory <- matrix(
    0.5, nrow = 4, ncol = 4,
    dimnames = list(reactions, units)
  )
  layer1 <- list(
    reaction_expression = multiome,
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
    "regcompass_corda_like_reaction_evidence_v1"
  )
})

test_that("CORDA-like classes retain HC and all module MC reactions", {
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
  expect_setequal(classes$biological, c("R1", "R2", "R3", "R4"))
})

test_that("NC support receives a larger weighted LP cost", {
  evidence <- data.frame(
    reaction_id = c("HC", "MC", "OT", "NC"),
    evidence_score = c(1, 0.8, 0.4, 0),
    stringsAsFactors = FALSE
  )
  classes <- RegCompassR:::.rc_corda_classify_reactions(
    parent_reactions = evidence$reaction_id,
    module_reactions = c("HC", "MC"),
    core_reactions = "HC",
    reaction_evidence = evidence,
    medium_confidence_threshold = 0.75,
    negative_confidence_threshold = 0.1,
    include_evidence_outside_modules = FALSE
  )
  costs <- RegCompassR:::.rc_corda_support_costs(
    evidence$reaction_id, classes,
    other_penalty = 1, negative_penalty = 10
  )
  expect_equal(costs[c("HC", "MC")], c(0, 0))
  expect_equal(costs[["OT"]], 1)
  expect_equal(costs[["NC"]], 10)
})

test_that("Layer 2 exposes optional CORDA-like completion without changing formals", {
  implementation <- paste(
    deparse(body(rc_regcompass_step_layer2)), collapse = "\n"
  )
  expect_match(implementation, "model_completion", fixed = TRUE)
  expect_match(implementation, "corda_like", fixed = TRUE)
  expect_match(
    implementation,
    ".rc_regcompass_step_layer2_completion_base",
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

test_that("CORDA-like cache falls back to the exact FASTCORE implementation", {
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
    "celltype_medium_corda_like_evidence_max",
    fixed = TRUE
  )
})
