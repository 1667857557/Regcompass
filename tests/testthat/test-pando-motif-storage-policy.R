test_that("RegCompass retains motif matches without exact hit positions", {
  active <- paste(
    deparse(body(.rc_run_condition_single_cell_grns)),
    collapse = "\n"
  )
  legacy <- paste(
    deparse(body(.rc_run_condition_single_cell_grns_legacy)),
    collapse = "\n"
  )
  multitask <- paste(
    deparse(body(.rc_run_celltype_multitask_grns_parameter_policy_core)),
    collapse = "\n"
  )

  for (workflow in list(active, legacy, multitask)) {
    expect_match(workflow, "store_motif_positions = FALSE", fixed = TRUE)
  }
})

test_that("motif position storage is workflow controlled", {
  expect_error(
    .rc_run_celltype_multitask_grns(
      object = NULL,
      gem = NULL,
      outdir = tempfile(),
      genome = NULL,
      pando_motif_args = list(store_motif_positions = TRUE)
    ),
    "store_motif_positions",
    fixed = TRUE
  )
})
