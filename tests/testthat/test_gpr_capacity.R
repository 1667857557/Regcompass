test_that("rc_reaction_capacity returns reaction by pool matrix", {
  gpr_table <- data.frame(
    reaction_id = c("R_HEX", "R_PFK"),
    gpr = c("HK1 or HK2", "PFKM and PFKL"),
    stringsAsFactors = FALSE
  )
  parsed <- rc_parse_gpr_table(gpr_table)
  gene_score <- matrix(
    c(0.8, 0.4, 0.7, 0.5, 0.2, 0.9, 0.6, 0.4),
    nrow = 4,
    dimnames = list(c("hk1", "hk2", "pfkm", "pfkl"), c("pool1", "pool2"))
  )

  out <- rc_reaction_capacity(parsed, gene_score, promiscuity_mode = "none", tau = 0.08)
  expect_equal(dim(out), c(2L, 2L))
  expect_identical(rownames(out), c("R_HEX", "R_PFK"))
  expect_identical(colnames(out), c("pool1", "pool2"))
  expect_equal(out["R_HEX", "pool1"], (0.8 + 0.4) / sqrt(2))

  sqrt_out <- rc_reaction_capacity(parsed, gene_score, promiscuity_mode = "sqrt", tau = 0.08)
  expect_equal(sqrt_out["R_HEX", "pool1"], out["R_HEX", "pool1"])
})

test_that("rc_reaction_capacity passes promiscuity_mode into weights", {
  gprs <- list(r1 = list("g1"), r2 = list("g1"))
  gene_score <- matrix(1, nrow = 1, ncol = 1, dimnames = list("g1", "p1"))
  none <- rc_reaction_capacity(gprs, gene_score, promiscuity_mode = "none")
  sqrt <- rc_reaction_capacity(gprs, gene_score, promiscuity_mode = "sqrt")
  linear <- rc_reaction_capacity(gprs, gene_score, promiscuity_mode = "linear")
  expect_equal(none["r1", "p1"], 1)
  expect_equal(sqrt["r1", "p1"], 1 / sqrt(2))
  expect_equal(linear["r1", "p1"], 1 / 2)
})

test_that("rc_run_layer1_capacity returns MVP Layer 1 toy outputs", {
  gpr_table <- data.frame(
    reaction_id = c("R_HEX", "R_PFK"),
    gpr = c("HK1 or HK2", "PFKM and PFKL"),
    stringsAsFactors = FALSE
  )
  expr <- matrix(
    c(10, 2, 5, 5, 4, 8, 6, 2),
    nrow = 4,
    dimnames = list(c("HK1", "HK2", "PFKM", "PFKL"), c("pool1", "pool2"))
  )
  detect <- matrix(
    c(1, 0.5, 0.8, 0.9, 1, 0.2, 0.7, 0.6),
    nrow = 4,
    dimnames = list(c("HK1", "HK2", "PFKM", "PFKL"), c("pool1", "pool2"))
  )

  out <- rc_run_layer1_capacity(gpr_table, expr, unit_detection = detect, promiscuity_mode = "sqrt", min_direct = 10, run_sensitivity = TRUE)
  expect_true(all(c("C_raw", "C_rel", "q95_diagnostics", "gpr_diagnostics", "reaction_confidence") %in% names(out)))
  expect_equal(dim(out$C_raw), c(2L, 2L))
  expect_true("mean_gpr_detection_rate" %in% colnames(out$reaction_confidence))

  pfk <- out$capacity_long[out$capacity_long$reaction_id == "R_PFK" & out$capacity_long$pool_id == "pool1", ]
  cap <- stats::setNames(pfk$C_raw, pfk$and_method)
  expect_lte(cap[["min"]], cap[["boltzmann_0.08"]])
  expect_lte(cap[["boltzmann_0.08"]], cap[["boltzmann_0.20"]])
  expect_lte(cap[["boltzmann_0.20"]], cap[["mean"]])
})

test_that("safe scale uses sigma-consistent MAD/IQR and clips z-scores", {
  x <- c(0, 0, 1, 1)
  expect_equal(rc_safe_scale(x, min_scale = 0.05), max(stats::mad(x, constant = 1.4826), stats::IQR(x) / 1.349, 0.05))
  X <- rbind(g1 = c(0, 100), g2 = c(1, 1))
  z <- rc_gene_zscore(X, z_clip = 2)
  expect_true(all(z <= 2 & z >= -2))
})

test_that("reaction confidence aggregates gene confidence by GPR AND/OR genes and pools", {
  gprs <- list(r1 = list(c("g1", "g2")), r2 = list("g3"))
  gene_conf <- matrix(c(0.2, 0.8, 0.9, 0.4, 0.6, 0.7), nrow = 3,
                      dimnames = list(c("g1", "g2", "g3"), c("p1", "p2")))
  out <- rc_reaction_confidence(gprs, gene_confidence = gene_conf, and_method = "min")
  expect_true(all(c("reaction_id", "pool_id", "reaction_confidence") %in% colnames(out)))
  expect_equal(out$reaction_confidence[out$reaction_id == "r1" & out$pool_id == "p1"], min(c(0.2, 0.8)))
})

test_that("reaction confidence falls back to detection when gene confidence has no GPR overlap", {
  gprs <- list(r1 = list(c("hk1", "pfkm")))
  gene_conf <- matrix(1, nrow = 1, ncol = 2, dimnames = list("scn5a", c("p1", "p2")))
  detect <- matrix(c(1, 0.5, 0.2, 0.8), nrow = 2, dimnames = list(c("hk1", "pfkm"), c("p1", "p2")))
  out <- rc_reaction_confidence(gprs, unit_detection = detect)
  expect_true(all(is.finite(out$reaction_confidence)))
  expect_true(all(out$detection_available))
})

test_that("metabolic GPR genes are extracted from the active GPR table", {
  gpr <- data.frame(reaction_id = c("r1", "r2"), gpr = c("HK1 or HK2", "PFKM and LDHA"))
  expect_setequal(rc_metabolic_gpr_genes(gpr), c("HK1", "HK2", "PFKM", "LDHA"))
})

test_that("rc_layer1_capacity alias matches run function", {
  expect_identical(rc_layer1_capacity, rc_run_layer1_capacity)
})

test_that("layer1 returns AND-method capacity long table and missing penalty", {
  gprs <- list(r1 = list(c("g1", "g2")))
  gene_score <- matrix(c(0.2, 0.8), nrow = 2, dimnames = list(c("g1", "g2"), "p1"))
  long <- rc_and_method_capacity_long(gprs, gene_score)
  expect_true(all(c("min", "boltzmann_0.08", "boltzmann_0.20", "mean") %in% long$and_method))
  conf <- rc_reaction_confidence(gprs, gene_confidence = matrix(1, nrow = 1, dimnames = list("g1", "p1")))
  expect_true(is.na(conf$reaction_confidence))
  expect_true(conf$reaction_unsupported_by_complete_gpr_flag)
})

test_that("rc_run_layer1_from_counts provides RNA-detection confidence source", {
  counts <- Matrix::Matrix(c(1, 0, 2, 3, 0, 4), nrow = 2, sparse = TRUE)
  rownames(counts) <- c("g1", "g2"); colnames(counts) <- paste0("c", 1:3)
  unit_map <- data.frame(pool_id = c("p1", "p1", "p2"), cell_id = colnames(counts), skipped = FALSE, sample_id = "s1", cell_type = "T")
  gpr <- data.frame(reaction_id = "r1", gpr = "g1 and g2")
  out <- rc_run_layer1_from_counts(gpr, counts, unit_map, bootstrap = FALSE)
  expect_equal(out$reaction_confidence_source, "gpr_aware_rna_detection")
  expect_true("capacity_long" %in% names(out))
})

test_that("layer1 ignores non-metabolic peak-gene links and keeps detection confidence", {
  counts <- Matrix::Matrix(c(1, 0, 2, 3, 0, 4), nrow = 2, sparse = TRUE)
  rownames(counts) <- c("HK1", "PFKM"); colnames(counts) <- paste0("c", 1:3)
  atac <- Matrix::Matrix(c(1, 1, 1), nrow = 1, sparse = TRUE)
  rownames(atac) <- "peak1"; colnames(atac) <- colnames(counts)
  unit_map <- data.frame(pool_id = paste0("p", 1:3), cell_id = colnames(counts), skipped = FALSE, sample_id = "s1", cell_type = "T")
  gpr <- data.frame(reaction_id = "r1", gpr = "HK1 and PFKM")
  links <- data.frame(peak_id = "peak1", gene = "SCN5A", weight = 1)
  out <- rc_run_layer1_from_counts(gpr, counts, unit_map, atac_counts = atac, peak_gene_links = links, bootstrap = FALSE)
  expect_equal(out$reaction_confidence_source, "gpr_aware_rna_detection")
  expect_true(all(is.finite(out$reaction_confidence$reaction_confidence)))
})

test_that("OR capacity supports bounded sensitivity methods", {
  x <- c(0.2, 0.5, 0.8)
  expect_equal(rc_or_capacity(x, method = "sum"), 1.5)
  expect_equal(rc_or_capacity(x, method = "max"), 0.8)
  expect_equal(rc_or_capacity(x, method = "prob_or"), 1 - prod(1 - x))
  expect_equal(rc_or_capacity(x, method = "sum_sqrtK"), sum(x) / sqrt(length(x)))
})

test_that("reaction confidence falls back per reaction when multiome overlap is partial", {
  gprs <- list(r_multi = list("g1"), r_detect = list("g2"))
  gene_conf <- matrix(c(0.8, 0.6), nrow = 1, dimnames = list("g1", c("p1", "p2")))
  detect <- matrix(c(0.2, 0.4, 0.7, 0.9), nrow = 2, dimnames = list(c("g1", "g2"), c("p1", "p2")))
  out <- rc_reaction_confidence(gprs, gene_confidence = gene_conf, unit_detection = detect)
  expect_equal(unique(out$confidence_source), "gpr_aware_gene_confidence")
  expect_true(all(is.finite(out$reaction_confidence[out$reaction_id == "r_multi"])))
  expect_true(all(is.na(out$reaction_confidence[out$reaction_id == "r_detect"])))
  expect_true(unique(out$reaction_unsupported_by_complete_gpr_flag[out$reaction_id == "r_detect"]))
})

test_that("Layer 1 names OR raw capacity according to OR method", {
  gpr_table <- data.frame(reaction_id = "r1", gpr = "g1 or g2", stringsAsFactors = FALSE)
  expr <- matrix(c(1, 3, 2, 4), nrow = 2, dimnames = list(c("g1", "g2"), c("p1", "p2")))
  out <- rc_run_layer1_capacity(gpr_table, expr, or_method = "max", run_sensitivity = FALSE)
  expect_equal(out$or_method_used, "max")
  expect_true("C_or_raw" %in% names(out))
  expect_false("C_iso_sum_raw" %in% names(out))
})

test_that("AND-method sensitivity keeps all-NA ranges as NA", {
  long <- data.frame(
    reaction_id = "r1",
    pool_id = "p1",
    and_method = c("min", "boltzmann_0.08", "boltzmann_0.20", "mean"),
    C_raw = NA_real_
  )
  out <- rc_and_method_sensitivity(long)
  expect_true(is.na(out$capacity_range))
  expect_false(out$tau_sensitive_flag)
})

test_that("missing AND subunits make capacity unavailable instead of reusing partial complex", {
  gprs <- list(r1 = list(c("g1", "g2")), r2 = list(c("g1", "g2"), "g3"))
  gene_score <- matrix(c(0.8, 0.4), nrow = 2, dimnames = list(c("g1", "g3"), "p1"))
  out <- rc_reaction_capacity(gprs, gene_score, promiscuity_mode = "none", or_method = "max")
  expect_true(is.na(out["r1", "p1"]))
  expect_equal(out["r2", "p1"], 0.4)
})

test_that("reaction confidence uses GPR bottleneck and isoenzyme semantics", {
  gprs <- list(r_and = list(c("g1", "g2")), r_or = list("g1", "g2"), r_missing = list(c("g1", "g3")))
  gene_conf <- matrix(c(0.2, 0.8), nrow = 2, dimnames = list(c("g1", "g2"), "p1"))
  out <- rc_reaction_confidence(gprs, gene_confidence = gene_conf, and_method = "min")
  expect_equal(out$reaction_confidence[out$reaction_id == "r_and"], 0.2)
  expect_equal(out$reaction_confidence[out$reaction_id == "r_or"], 0.8)
  expect_true(is.na(out$reaction_confidence[out$reaction_id == "r_missing"]))
  expect_true(out$reaction_unsupported_by_complete_gpr_flag[out$reaction_id == "r_missing"])
})

test_that("GPR-aware confidence does not penalize supported OR isoenzymes", {
  gprs <- list(R1 = list("A", "B", "C"))
  gene_conf <- matrix(c(0.9, 0.0, 0.0), nrow = 3, dimnames = list(c("A", "B", "C"), "p1"))
  aware <- rc_reaction_confidence_gpr_aware(gprs, gene_confidence = gene_conf, or_method = "max")
  expect_equal(aware$reaction_confidence, 0.9)
  expect_false(aware$no_complete_gpr_group_flag)
})

test_that("GPR-aware confidence limits AND complexes by low subunits", {
  gprs <- list(R2 = list(c("A", "B", "C")))
  gene_conf <- matrix(c(0.9, 0.8, 0.1), nrow = 3, dimnames = list(c("A", "B", "C"), "p1"))
  aware_min <- rc_reaction_confidence_gpr_aware(gprs, gene_confidence = gene_conf, and_method = "min")
  aware_soft <- rc_reaction_confidence_gpr_aware(gprs, gene_confidence = gene_conf, and_method = "softmin", tau_conf = 0.20)
  expect_equal(aware_min$reaction_confidence, 0.1)
  expect_lt(aware_soft$reaction_confidence, 0.2)
})

test_that("GPR-aware confidence uses complete alternative AND groups", {
  gprs <- list(R3 = list(c("A", "B"), "C"))
  gene_conf <- matrix(0.85, nrow = 1, dimnames = list("C", "p1"))
  aware <- rc_reaction_confidence_gpr_aware(gprs, gene_confidence = gene_conf)
  expect_equal(aware$reaction_confidence, 0.85)
  expect_false(aware$no_complete_gpr_group_flag)
  expect_equal(aware$n_and_groups_complete, 1L)
})

test_that("GPR-aware confidence keeps fully missing reactions as NA", {
  gprs <- list(R4 = list(c("X", "Y")))
  gene_conf <- matrix(0.5, nrow = 1, dimnames = list("A", "p1"))
  aware <- rc_reaction_confidence_gpr_aware(gprs, gene_confidence = gene_conf)
  expect_true(is.na(aware$reaction_confidence))
  expect_true(aware$no_complete_gpr_group_flag)
})

test_that("RNA-only GPR-aware confidence does not median-penalize OR isoenzymes", {
  gprs <- list(R1 = list(c("A"), c("B"), c("C")))
  det <- matrix(c(0.9, 0, 0), nrow = 3, dimnames = list(c("a", "b", "c"), "pool1"))
  out <- rc_reaction_confidence(gprs, unit_detection = det, method = "gpr_aware")
  expect_equal(out$reaction_confidence, 0.9)
  expect_equal(out$confidence_source, "gpr_aware_rna_detection")
})

test_that("GPR-aware confidence distinguishes incomplete alternatives from unsupported reactions", {
  gprs <- list(R3 = list(c("A", "B"), c("C")))
  ev <- matrix(c(NA_real_, NA_real_, 0.85), nrow = 3, dimnames = list(c("a", "b", "c"), "pool1"))
  out <- rc_reaction_confidence_gpr_aware(gprs, ev)
  expect_false(out$no_complete_gpr_group_flag)
  expect_false(out$reaction_unsupported_by_complete_gpr_flag)
  expect_equal(out$reaction_confidence, 0.85)
})

test_that("rc_run_layer1_capacity uses GPR-aware RNA-only detection confidence", {
  gpr_table <- data.frame(reaction_id = "R1", gpr = "A or B or C", stringsAsFactors = FALSE)
  expr <- matrix(c(1, 1, 1), nrow = 3, dimnames = list(c("A", "B", "C"), "pool1"))
  det <- matrix(c(0.9, 0, 0), nrow = 3, dimnames = list(c("A", "B", "C"), "pool1"))
  out <- rc_run_layer1_capacity(gpr_table, expr, unit_detection = det, bootstrap = FALSE)
  expect_equal(out$reaction_confidence$reaction_confidence, 0.9)
  expect_equal(out$reaction_confidence$confidence_source, "gpr_aware_rna_detection")
})
