test_that("edge ordering precedes predictor extraction", {
  body_text <- paste(
    deparse(body(.rc_fit_multitask_target_direct)),
    collapse = "\n"
  )
  order_position <- regexpr(".rc_order_target_edges", body_text, fixed = TRUE)[1]
  rna_position <- regexpr("rna[edges$tf_feature_id", body_text, fixed = TRUE)[1]
  atac_position <- regexpr("atac[edges$atac_feature_id", body_text, fixed = TRUE)[1]
  expect_gt(order_position, 0)
  expect_gt(rna_position, order_position)
  expect_gt(atac_position, order_position)
})
