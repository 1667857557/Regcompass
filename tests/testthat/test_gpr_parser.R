test_that("rc_parse_gpr_simple supports v0.3 grammar", {
  expect_identical(rc_parse_gpr_simple("(HK1 and HK2) or HK3"), list(c("hk1", "hk2"), "hk3"))
  expect_identical(rc_parse_gpr_simple("PFKM and PFKL"), list(c("pfkm", "pfkl")))
  expect_identical(rc_parse_gpr_simple("LDHA or LDHB"), list("ldha", "ldhb"))
  expect_identical(rc_parse_gpr_simple("GAPDH"), list("gapdh"))
})

test_that("rc_promiscuity_weight supports sensitivity modes", {
  gprs <- list(
    r1 = list(c("g1", "g2")),
    r2 = list("g1"),
    r3 = list("g3")
  )
  expect_equal(rc_promiscuity_weight(gprs, "none")["g1"], c(g1 = 1))
  expect_equal(rc_promiscuity_weight(gprs, "sqrt")["g1"], c(g1 = 1 / sqrt(2)))
  expect_equal(rc_promiscuity_weight(gprs, "linear")["g1"], c(g1 = 1 / 2))
})

test_that("rc_parse_gpr_table supports long-table GPR input", {
  tab <- data.frame(reaction_id = c("r1", "r1", "r1"), and_group_id = c(1, 1, 2), gene = c("G1", "G2", "G3"))
  expect_identical(rc_parse_gpr_table(tab), list(r1 = list(c("g1", "g2"), "g3")))
})

test_that("rc_parse_gpr_simple rejects complex nested GPR input", {
  expect_error(rc_parse_gpr_simple("g1 and (g2 or g3)"), "Complex")
})

test_that("rc_parse_gpr_simple supports flat parenthesized OR-of-AND", {
  expect_identical(rc_parse_gpr_simple("(g1 and g2) or (g3 and g4)"), list(c("g1", "g2"), c("g3", "g4")))
})
