test_that("discordant shared-TF edges do not merge biological components", {
  input <- data.frame(
    sample_id = c("s1", "s1"),
    tf = c("TF1", "TF1"),
    target = c("G1", "G2"),
    estimate = c(1, -2),
    stringsAsFactors = FALSE
  )
  projected <- rc_project_metabolic_grn(
    input,
    metabolic_genes = c("G1", "G2"),
    top_k = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0
  )
  expect_equal(projected$edges$regulatory_relation, "discordant")
  expect_false(projected$edges$used_for_component)
  expect_equal(length(unique(projected$nodes$module_id)), 2L)
})

test_that("concordant shared-TF edges can define one component", {
  input <- data.frame(
    sample_id = c("s1", "s1"),
    tf = c("TF1", "TF1"),
    target = c("G1", "G2"),
    estimate = c(1, 2),
    stringsAsFactors = FALSE
  )
  projected <- rc_project_metabolic_grn(
    input,
    metabolic_genes = c("G1", "G2"),
    top_k = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0
  )
  expect_equal(projected$edges$regulatory_relation, "concordant")
  expect_true(projected$edges$used_for_component)
  expect_equal(length(unique(projected$nodes$module_id)), 1L)
})
