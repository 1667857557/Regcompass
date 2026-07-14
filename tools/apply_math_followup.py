from pathlib import Path


def replace_once(path, old, new):
    p = Path(path)
    text = p.read_text()
    if old not in text:
        raise RuntimeError(f"Expected block not found in {path}")
    p.write_text(text.replace(old, new, 1))


replace_once(
    "R/stats.R",
    '''.rc_aggregate_microcompass_samples <- function(score, meta, sample_col,
                                                condition_col, covariates = NULL) {
  score <- as.matrix(score)
  fields <- unique(c(sample_col, condition_col, covariates))
  missing <- setdiff(fields, colnames(meta))
  if (length(missing)) {
    stop("Metadata is missing sample-level fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  samples <- unique(as.character(meta[[sample_col]]))
  if (anyNA(samples) || any(!nzchar(samples))) {
    stop("Biological sample IDs must be non-missing and non-empty.", call. = FALSE)
  }
  sample_meta <- lapply(samples, function(sample) {
    rows <- which(as.character(meta[[sample_col]]) == sample)
    values <- lapply(fields, function(field) unique(as.character(meta[[field]][rows])))
    names(values) <- fields
    bad <- vapply(values, function(value) length(value) != 1L || is.na(value) || !nzchar(value), logical(1))
    if (any(bad)) {
      stop("Sample `", sample, "` has inconsistent metadata: ",
           paste(fields[bad], collapse = ", "), call. = FALSE)
    }
    as.data.frame(c(list(sample_id_internal = sample), values), stringsAsFactors = FALSE)
  })
  sample_meta <- do.call(rbind, sample_meta)
  rownames(sample_meta) <- sample_meta$sample_id_internal
  sample_meta$sample_id_internal <- NULL
  aggregated <- do.call(cbind, lapply(samples, function(sample) {
    columns <- which(as.character(meta[[sample_col]]) == sample)
    matrixStats::rowMedians(score[, columns, drop = FALSE], na.rm = TRUE)
  }))
  rownames(aggregated) <- rownames(score)
  colnames(aggregated) <- samples
  list(score = aggregated, meta = sample_meta)
}
''',
    '''.rc_aggregate_microcompass_samples <- function(score, meta, sample_col,
                                                condition_col, covariates = NULL) {
  score <- as.matrix(score)
  fields <- unique(c(sample_col, condition_col, covariates))
  missing <- setdiff(fields, colnames(meta))
  if (length(missing)) {
    stop("Metadata is missing sample-level fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  samples <- unique(as.character(meta[[sample_col]]))
  if (anyNA(samples) || any(!nzchar(samples))) {
    stop("Biological sample IDs must be non-missing and non-empty.", call. = FALSE)
  }
  sample_meta <- lapply(samples, function(sample) {
    rows <- which(as.character(meta[[sample_col]]) == sample)
    bad <- vapply(fields, function(field) {
      values <- meta[[field]][rows]
      values <- values[!is.na(values)]
      length(unique(values)) != 1L
    }, logical(1))
    if (any(bad)) {
      stop("Sample `", sample, "` has inconsistent metadata: ",
           paste(fields[bad], collapse = ", "), call. = FALSE)
    }
    meta[rows[[1L]], fields, drop = FALSE]
  })
  sample_meta <- do.call(rbind, sample_meta)
  rownames(sample_meta) <- samples
  aggregated <- do.call(cbind, lapply(samples, function(sample) {
    columns <- which(as.character(meta[[sample_col]]) == sample)
    matrixStats::rowMedians(score[, columns, drop = FALSE], na.rm = TRUE)
  }))
  rownames(aggregated) <- rownames(score)
  colnames(aggregated) <- samples
  list(score = aggregated, meta = sample_meta)
}
'''
)

replace_once(
    "R/microcompass.R",
    '''  sample_conditions <- .rc_sample_condition_map(
    matrices$unit_meta, sample_col, condition_col
  )

  if (identical(mode, "full_gem")) {''',
    '''  sample_conditions <- .rc_sample_condition_map(
    matrices$unit_meta, sample_col, condition_col
  )
  direction_diagnostics <- NULL

  if (identical(mode, "full_gem")) {'''
)

replace_once(
    "R/microcompass.R",
    '''    directions <- rc_prepare_directional_targets(
      gem,
      target_reactions,
      target_direction
    )
    if (!nrow(directions)) {
      stop(
        "No target reaction directions are allowed by the GEM bounds.",
        call. = FALSE
      )
    }
''',
    '''    directions <- rc_prepare_directional_targets(
      gem,
      target_reactions,
      target_direction
    )
    direction_diagnostics <- directions
    directions <- directions[
      directions$target_direction %in% c("forward", "reverse"),
      , drop = FALSE
    ]
    if (!nrow(directions)) {
      stop(
        "No target reaction directions are allowed by the GEM bounds.",
        call. = FALSE
      )
    }
'''
)

replace_once(
    "R/microcompass.R",
    '''    target_direction = directions,
    medium_scenarios = medium_scenarios,''',
    '''    target_direction = directions,
    direction_diagnostics = direction_diagnostics,
    medium_scenarios = medium_scenarios,'''
)

test_file = Path("tests/testthat/test_math_biology_audit_v13.R")
test_text = test_file.read_text()
addition = r'''

test_that("sample aggregation preserves continuous covariate classes", {
  score <- matrix(1:8, nrow = 1,
                  dimnames = list("R::forward::medium=base", paste0("u", 1:8)))
  meta <- data.frame(
    unit_id = paste0("u", 1:8),
    sample_id = rep(c("S1", "S2", "S3", "S4"), each = 2),
    condition = rep(c("A", "A", "B", "B"), each = 2),
    age = rep(c(30, 40, 50, 60), each = 2),
    stringsAsFactors = FALSE
  )
  out <- .rc_aggregate_microcompass_samples(
    score, meta, "sample_id", "condition", covariates = "age"
  )
  expect_true(is.numeric(out$meta$age))
  expect_equal(out$meta$age, c(30, 40, 50, 60))
})

test_that("blocked full-GEM directions are retained only as diagnostics", {
  S <- matrix(0, nrow = 1, ncol = 2,
              dimnames = list("m", c("blocked", "open")))
  gem <- rc_make_gem(S, lb = c(blocked = 0, open = 0),
                     ub = c(blocked = 0, open = 10))
  directions <- rc_prepare_directional_targets(
    gem, c("blocked", "open"), target_direction = "both"
  )
  allowed <- directions[directions$target_direction %in% c("forward", "reverse"), , drop = FALSE]
  expect_equal(directions$target_direction[directions$reaction_id == "blocked"], "none")
  expect_false("blocked" %in% allowed$reaction_id)
  expect_equal(allowed$target_direction, "forward")
})
'''
if 'sample aggregation preserves continuous covariate classes' not in test_text:
    test_file.write_text(test_text.rstrip() + addition + "\n")
