test_that("CORDA options preserve FASTCORE and match paper defaults", {
  defaults <- RegCompassR:::.rc_layer2_corda_options(list())
  expect_identical(defaults$model_completion, "fastcore")
  corda <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda"
  ))
  expect_identical(corda$model_completion, "corda")
  expect_equal(corda$gamma, 1e5)
  expect_equal(corda$kappa, 1e-2)
  expect_equal(corda$epsilon, 1)
  expect_equal(corda$n, 5L)
  expect_equal(corda$p, 2L)
  alias <- RegCompassR:::.rc_layer2_corda_options(list(
    model_completion = "corda_like"
  ))
  expect_identical(alias$model_completion, "corda")
  expect_error(
    RegCompassR:::.rc_layer2_corda_options(list(
      model_completion = "corda",
      corda_negative_penalty = 10
    )),
    "Obsolete weighted-FASTCORE"
  )
})

test_that("CORDA evidence is calculated within each cell type", {
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

test_that("confidence mapping leaves MC flexible", {
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
  expect_setequal(classes$mc, c("R2", "R3", "R4"))
  expect_identical(classes$nc, "R6")
  expect_setequal(classes$ot, c("R5", "R7"))
  expect_false("biological" %in% names(classes))
  expect_true(all(classes$confidence[classes$mc] != "RE"))
})

test_that("direction splitting preserves bounds and direction", {
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
  expect_setequal(
    split$direction_table$variable_id,
    c("IRR::forward", "REV::forward", "REV::reverse")
  )
  expect_equal(split$ub[["IRR::forward"]], 10)
  expect_equal(split$ub[["REV::forward"]], 7)
  expect_equal(split$ub[["REV::reverse"]], 5)
  expect_equal(
    as.numeric(split$S[, "REV::reverse"]),
    -as.numeric(S[, "REV"])
  )
})

test_that("Layer 2 exposes CORDA without changing public formals", {
  implementation <- paste(
    deparse(body(rc_regcompass_step_layer2)), collapse = "\n"
  )
  expect_match(implementation, "corda", fixed = TRUE)
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

test_that("CORDA cache retains exact FASTCORE fallback", {
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
    "celltype_medium_original_corda",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "models_serial_target_direction_x_replicate_inner_parallel",
    fixed = TRUE
  )
})
