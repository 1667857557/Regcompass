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
