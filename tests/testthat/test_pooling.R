test_that("rc_make_pools never mixes samples or cell types", {
  meta <- data.frame(
    sample_id = rep(c("s1", "s2"), each = 60),
    condition = rep(c("case", "control"), each = 60),
    cell_type = rep(rep(c("T", "B"), each = 30), 2),
    seurat_clusters = rep(c("0", "1"), each = 15, times = 4),
    row.names = paste0("cell", seq_len(120)),
    stringsAsFactors = FALSE
  )

  pools <- rc_make_pools(
    meta,
    condition_col = "condition",
    state_col = "seurat_clusters",
    target_size = 20,
    min_size = 10,
    seed = 42
  )

  expect_true(all(c("pool_id", "cell_id", "low_power_pool", "pool_size", "group_key") %in% colnames(pools)))
  expect_identical(sort(pools$cell_id), sort(rownames(meta)))

  per_pool <- split(pools, pools$pool_id)
  expect_true(all(vapply(per_pool, function(x) length(unique(x$sample_id)) == 1L, logical(1))))
  expect_true(all(vapply(per_pool, function(x) length(unique(x$cell_type)) == 1L, logical(1))))
  expect_true(all(vapply(per_pool, function(x) length(unique(x$condition)) == 1L, logical(1))))
  expect_true(all(vapply(per_pool, function(x) length(unique(x$seurat_clusters)) == 1L, logical(1))))
})

test_that("rc_make_pools falls back to sample by cell type when no state is supplied", {
  meta <- data.frame(
    sample_id = rep(c("s1", "s2"), each = 35),
    cell_type = rep(c("T", "B"), each = 35),
    row.names = paste0("cell", seq_len(70)),
    stringsAsFactors = FALSE
  )

  pools <- rc_make_pools(meta, target_size = 20, min_size = 10, seed = 1)
  expect_false("seurat_clusters" %in% colnames(pools))
  expect_equal(length(unique(pools$group_key)), 2L)
  expect_true(all(vapply(split(pools, pools$pool_id), function(x) length(unique(x$sample_id)) == 1L, logical(1))))
})

test_that("rc_make_pools flags low-power strata and validates columns", {
  meta <- data.frame(
    sample_id = "s1",
    cell_type = "rare",
    row.names = paste0("cell", seq_len(5)),
    stringsAsFactors = FALSE
  )

  pools <- rc_make_pools(meta, target_size = 80, min_size = 30, seed = 1)
  expect_equal(length(unique(pools$pool_id)), 1L)
  expect_true(all(pools$low_power_pool))
  expect_equal(unique(pools$pool_size), 5L)

  expect_error(rc_make_pools(meta, condition_col = "condition"), "Missing metadata columns: condition")
})

test_that("pool seed replicates record seed metadata", {
  meta <- data.frame(sample_id = "s1", cell_type = "T", row.names = paste0("c", 1:40))
  pools <- rc_make_pool_seed_replicates(meta, seeds = c(1, 2), target_size = 20, min_group_size = 10, min_pool_size = 5)
  expect_equal(sort(unique(pools$pool_seed_replicate)), c(1, 2))
  expect_true("pool_seed" %in% colnames(pools))
  expect_equal(anyDuplicated(stats::na.omit(pools$pool_id)), 0L)
})

test_that("rc_make_pools uses harmony or pca embeddings when requested", {
  skip_if_not_installed("SeuratObject")

  meta <- data.frame(
    sample_id = "s1",
    cell_type = "T",
    row.names = paste0("cell", seq_len(60)),
    stringsAsFactors = FALSE
  )
  counts <- matrix(0, nrow = 3, ncol = 60, dimnames = list(paste0("g", 1:3), rownames(meta)))
  seu <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  emb <- matrix(rnorm(60 * 3), nrow = 60, dimnames = list(rownames(meta), paste0("PC_", 1:3)))
  seu[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "PC_", assay = "RNA")

  pools <- rc_make_pools(
    meta = meta,
    seu = seu,
    target_size = 30,
    min_group_size = 10,
    min_pool_size = 5,
    pooling_method = "auto",
    dims = 1:2,
    seed = 1
  )

  expect_equal(unique(pools$pooling_method), "embedding")
  expect_equal(unique(pools$pool_reduction), "pca")
  expect_equal(unique(pools$pool_dims), "1,2")
  expect_true(all(vapply(split(pools, pools$pool_id), function(x) length(unique(x$sample_id)) == 1L, logical(1))))
})
