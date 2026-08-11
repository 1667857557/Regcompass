test_that("condition target worker budgets evenly exhaust the cap", {
  expect_identical(
    .rc_allocate_condition_target_workers(60L, 8L),
    c(8L, 8L, 8L, 8L, 7L, 7L, 7L, 7L)
  )
  expect_equal(sum(.rc_allocate_condition_target_workers(60L, 8L)), 60L)
  expect_identical(
    .rc_allocate_condition_target_workers(10L, 3L),
    c(4L, 3L, 3L)
  )
})

test_that("more condition GRNs than workers stay one-worker tasks", {
  budget <- .rc_allocate_condition_target_workers(4L, 8L)
  expect_identical(budget, rep(1L, 8L))
  expect_equal(sum(head(budget, 4L)), 4L)
})

test_that("nested condition target backend is bounded", {
  skip_if_not_installed("BiocParallel")
  expect_null(.rc_condition_nested_target_bpparam(1L))
  param <- .rc_condition_nested_target_bpparam(3L)
  expect_s4_class(param, "SnowParam")
  expect_equal(BiocParallel::bpnworkers(param), 3L)
  expect_equal(attr(param, "regcompass_worker_limit"), 3L)
})
