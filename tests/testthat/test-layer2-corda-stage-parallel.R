.make_corda_stage_parallel_case <- function() {
  metabolites <- c("A1", "B1", "A2", "B2", "C1", "C2")
  reactions <- c("M1", "H1", "N1", "M2", "H2", "N2", "M3", "N3", "M4", "N4")
  S <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(length(metabolites), length(reactions)), dimnames = list(metabolites, reactions), giveCsparse = TRUE)
  S["A1", "M1"] <- 1; S["A1", "H1"] <- -1; S["B1", "H1"] <- 1; S["B1", "N1"] <- -1
  S["A2", "M2"] <- 1; S["A2", "H2"] <- -1; S["B2", "H2"] <- 1; S["B2", "N2"] <- -1
  S["C1", "M3"] <- 1; S["C1", "N3"] <- -1; S["C2", "M4"] <- 1; S["C2", "N4"] <- -1
  split <- RegCompassR:::.rc_corda2_split_original(list(S = S, lb = stats::setNames(rep(0, length(reactions)), reactions), ub = stats::setNames(rep(1000, length(reactions)), reactions)))
  initial <- stats::setNames(rep("OT", length(reactions)), reactions)
  initial[c("H1", "H2")] <- "HC"; initial[c("M1", "M2", "M3", "M4")] <- "MC_module"; initial[c("N1", "N2", "N3", "N4")] <- "NC"
  classes <- list(hc = c("H1", "H2"), mc_module = c("M1", "M2", "M3", "M4"), mc_evidence = character(), mc = c("M1", "M2", "M3", "M4"), nc = c("N1", "N2", "N3", "N4"), ot = character(), confidence = initial, initial_confidence = initial)
  list(split = split, classes = classes, options = RegCompassR:::.rc_layer2_corda_options(list(model_completion = "corda2")))
}

test_that("stage-parallel CORDA2 preserves exact serial mathematical state", {
  skip_if_not_installed("BiocParallel")
  case <- .make_corda_stage_parallel_case(); old_options <- options(RegCompassR.progress = FALSE); on.exit(options(old_options), add = TRUE)
  previous <- RegCompassR:::.rc_layer2_enter_parallel_context(FALSE, FALSE)
  serial <- RegCompassR:::.rc_corda_build_three_stage_core(split = case$split, classes = case$classes, options = case$options, solver = "highs", time_limit = 30)
  RegCompassR:::.rc_layer2_restore_parallel_context(previous)
  param <- BiocParallel::SnowParam(workers = 2L, type = "SOCK", progressbar = FALSE)
  previous <- RegCompassR:::.rc_layer2_enter_parallel_context(TRUE, param)
  on.exit({ RegCompassR:::.rc_layer2_restore_parallel_context(previous); try(BiocParallel::bpstop(param), silent = TRUE) }, add = TRUE)
  parallel <- RegCompassR:::.rc_corda_build_three_stage_core(split = case$split, classes = case$classes, options = case$options, solver = "highs", time_limit = 30)
  expect_setequal(parallel$included, serial$included); expect_setequal(parallel$included_directional_variables, serial$included_directional_variables)
  expect_equal(parallel$HCtoMC, serial$HCtoMC); expect_equal(parallel$HCtoNC, serial$HCtoNC); expect_equal(parallel$MCtoNC, serial$MCtoNC)
  expect_setequal(parallel$stage1_associated, serial$stage1_associated); expect_setequal(parallel$stage2_promoted_nc, serial$stage2_promoted_nc); expect_setequal(parallel$stage2_promoted_mc, serial$stage2_promoted_mc); expect_setequal(parallel$stage3_associated_ot, serial$stage3_associated_ot)
  expect_identical(serial$parallel_execution_policy, "serial_original_persistent_engine")
  expect_identical(parallel$parallel_execution_policy, "stage_barrier_parallel_targets_deterministic_ordered_reduce")
  expect_false(BiocParallel::bpisup(param))
})

test_that("canonical CORDA2 core keeps serial route when parallel is disabled", {
  case <- .make_corda_stage_parallel_case(); previous <- RegCompassR:::.rc_layer2_enter_parallel_context(FALSE, FALSE); on.exit(RegCompassR:::.rc_layer2_restore_parallel_context(previous), add = TRUE)
  result <- RegCompassR:::.rc_corda_build_three_stage_core(split = case$split, classes = case$classes, options = case$options, solver = "highs", time_limit = 30)
  expect_identical(result$parallel_execution_policy, "serial_original_persistent_engine")
})
