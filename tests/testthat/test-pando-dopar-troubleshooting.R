test_that("Pando dependency pin and dopar troubleshooting stay aligned", {
  description <- read.dcf(test_path("..", "..", "DESCRIPTION"))
  readme <- paste(
    readLines(test_path("..", "..", "README.md"), warn = FALSE),
    collapse = "\n"
  )
  pin <- "1667857557/Pando_regcompass@6f42c8143bec6610b001e714a51627337f6d9ba9"

  expect_match(description[[1L, "Remotes"]], pin, fixed = TRUE)
  expect_match(readme, pin, fixed = TRUE)
  expect_match(readme, "could not find function \"%dopar%\"", fixed = TRUE)
  expect_match(readme, "table(A$dataset, useNA = \"ifany\")", fixed = TRUE)
})
