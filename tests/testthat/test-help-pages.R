test_that("every exported function has an installed help topic", {
  exports <- getNamespaceExports("RegCompassR")
  for (name in exports) {
    topic <- utils::help(name, package = "RegCompassR")
    expect_length(topic, 1L, info = paste("missing help for", name))
  }
})
