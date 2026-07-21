test_that("offline Pando installation without remote metadata is allowed", {
  offline_description <- list(
    Package = "Pando",
    Version = "1.0.0"
  )

  expect_warning(
    result <- .rc_validate_pando_repository(
      description = offline_description,
      installed_version = "1.0.0"
    ),
    "offline or local source-package installation"
  )

  expect_identical(result$version, "1.0.0")
  expect_true(is.na(result$remote_username))
  expect_true(is.na(result$remote_repo))
  expect_false(result$repository_verified)
  expect_identical(
    result$installation_source,
    "local_or_offline_source_unverified"
  )
})

test_that("empty Pando remote fields are treated as offline metadata", {
  offline_description <- list(
    RemoteUsername = " ",
    RemoteRepo = ""
  )

  expect_warning(
    result <- .rc_validate_pando_repository(
      description = offline_description,
      installed_version = "1.0.0"
    ),
    "remote metadata are unavailable"
  )
  expect_false(result$repository_verified)
})

test_that("partial or conflicting Pando remote metadata still fail", {
  expect_error(
    .rc_validate_pando_repository(
      description = list(RemoteRepo = "Pando_regcompass"),
      installed_version = "1.0.0"
    ),
    "remote username mismatch"
  )

  expect_error(
    .rc_validate_pando_repository(
      description = list(
        RemoteUsername = "1667857557",
        RemoteRepo = "Pando"
      ),
      installed_version = "1.0.0"
    ),
    "remote repository mismatch"
  )
})

test_that("verified Pando remote metadata remain strict", {
  result <- .rc_validate_pando_repository(
    description = list(
      RemoteUsername = "1667857557",
      RemoteRepo = "Pando_regcompass",
      RemoteRef = "main",
      RemoteSha = "abc123"
    ),
    installed_version = "1.0.0"
  )

  expect_true(result$repository_verified)
  expect_identical(result$installation_source, "github_remote_verified")
  expect_identical(result$remote_ref, "main")
  expect_identical(result$remote_sha, "abc123")
})
