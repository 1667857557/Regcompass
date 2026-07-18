from pathlib import Path

def replace_once(path: str, old: str, new: str) -> None:
    file_path = Path(path)
    text = file_path.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one match, found {count}")
    file_path.write_text(text.replace(old, new, 1), encoding="utf-8")

# Bump bug-fix version.
replace_once(
    "DESCRIPTION",
    "Version: 1.5.0\n",
    "Version: 1.5.1\n",
)

# Invalidate incompatible species GEM caches and rebuild them automatically.
replace_once(
    "R/humangem.R",
    '''}

#' Prepare a species-specific genome-scale metabolic model
''',
    '''}

.rc_load_compatible_species_gem <- function(save_rds, spec) {
  if (!file.exists(save_rds)) return(NULL)
  cached <- tryCatch(
    rc_read_gem(save_rds),
    error = function(error) error
  )
  reason <- NULL
  if (inherits(cached, "error")) {
    reason <- conditionMessage(cached)
  } else {
    reason <- tryCatch(
      {
        rc_validate_species_gem(cached, spec$species)
        recorded_source <- as.character(cached$model_info$source %||% "")
        recorded_version <- as.character(cached$model_info$version %||% "")
        if (!identical(recorded_source, spec$source)) {
          stop(
            "cached model source is `", recorded_source,
            "` instead of `", spec$source, "`",
            call. = FALSE
          )
        }
        if (!identical(recorded_version, spec$version)) {
          stop(
            "cached model version is `", recorded_version,
            "` instead of `", spec$version, "`",
            call. = FALSE
          )
        }
        NULL
      },
      error = function(error) conditionMessage(error)
    )
  }
  if (!is.null(reason)) {
    warning(
      "Removing incompatible cached ", spec$repository_name,
      " model at `", save_rds, "`: ", reason,
      call. = FALSE
    )
    unlink(save_rds, force = TRUE)
    return(NULL)
  }
  cached
}

#' Prepare a species-specific genome-scale metabolic model
''',
)

old_prepare_block = '''  if (isTRUE(force_download) && file.exists(save_rds)) {
    unlink(save_rds, force = TRUE)
  }
  if (!file.exists(save_rds)) {
    ref <- if (identical(spec$version, "latest")) {
      "main"
    } else {
      paste0("v", spec$version)
    }
    tmp <- tempfile(paste0(spec$repository_name, "-"))
    prepared <- rc_download_species_gem(
      species = species,
      destdir = tmp,
      ref = ref,
      ref_type = "auto",
      gene_format = spec$gene_format,
      overwrite = TRUE,
      quiet = TRUE
    )
    repo_dir <- attr(prepared, "repo_dir")
    model_yml <- file.path(repo_dir, "model", spec$model_file)
    checksum <- unname(tools::md5sum(model_yml)[[1L]])
    gem <- rc_convert_yaml_to_regcompass(
      model_yml = model_yml,
      species = species,
      version = spec$version,
      commit = ref,
      checksum = checksum
    )
    gem <- rc_enrich_humangem_metadata(
      gem,
      reactions_tsv = prepared$reactions,
      model_yml = model_yml
    )
    gem <- rc_annotate_reaction_roles(gem, overwrite_existing = TRUE)
    gem$gpr_table <- prepared$gpr_table
    gem$metabolic_genes <- prepared$metabolic_genes
    gem$reaction_rules <- prepared$reaction_rules
    gem$genes <- prepared$genes
    gem$reactions <- prepared$reactions
    gem$model_info$gene_format <- spec$gene_format
    gem$model_info$archive <- attr(prepared, "archive")
    gem$model_info$archive_url <- attr(prepared, "archive_url") %||%
      NA_character_
    gem$model_info$annotation_schema <- "regcompass_species_gem_v1"
    gem$model_info$citation <- spec$citation
    gem$model_info$citation_doi <- spec$citation_doi
    rc_validate_species_gem(gem, species)
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(gem, save_rds)
  }
  gem <- rc_read_gem(save_rds)
'''
new_prepare_block = '''  if (isTRUE(force_download) && file.exists(save_rds)) {
    unlink(save_rds, force = TRUE)
  }
  cached <- if (isTRUE(force_download)) {
    NULL
  } else {
    .rc_load_compatible_species_gem(save_rds, spec)
  }
  if (!is.null(cached)) return(cached)

  ref <- if (identical(spec$version, "latest")) {
    "main"
  } else {
    paste0("v", spec$version)
  }
  tmp <- tempfile(paste0(spec$repository_name, "-"))
  prepared <- rc_download_species_gem(
    species = species,
    destdir = tmp,
    ref = ref,
    ref_type = "auto",
    gene_format = spec$gene_format,
    overwrite = TRUE,
    quiet = TRUE
  )
  repo_dir <- attr(prepared, "repo_dir")
  model_yml <- file.path(repo_dir, "model", spec$model_file)
  checksum <- unname(tools::md5sum(model_yml)[[1L]])
  gem <- rc_convert_yaml_to_regcompass(
    model_yml = model_yml,
    species = species,
    version = spec$version,
    commit = ref,
    checksum = checksum
  )
  gem <- rc_enrich_humangem_metadata(
    gem,
    reactions_tsv = prepared$reactions,
    model_yml = model_yml
  )
  gem <- rc_annotate_reaction_roles(gem, overwrite_existing = TRUE)
  gem$gpr_table <- prepared$gpr_table
  gem$metabolic_genes <- prepared$metabolic_genes
  gem$reaction_rules <- prepared$reaction_rules
  gem$genes <- prepared$genes
  gem$reactions <- prepared$reactions
  gem$model_info$gene_format <- spec$gene_format
  gem$model_info$archive <- attr(prepared, "archive")
  gem$model_info$archive_url <- attr(prepared, "archive_url") %||%
    NA_character_
  gem$model_info$annotation_schema <- "regcompass_species_gem_v1"
  gem$model_info$citation <- spec$citation
  gem$model_info$citation_doi <- spec$citation_doi
  rc_validate_species_gem(gem, species)
  dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(gem, save_rds)
  gem <- rc_read_gem(save_rds)
'''
replace_once("R/humangem.R", old_prepare_block, new_prepare_block)

# Key full-GEM cache files by the actual GEM structure and stable provenance.
replace_once(
    "R/full_gem.R",
    '''#' Cache one complete full GEM per medium scenario
rc_build_full_gem_cache <- function(gem, dirs, medium_scenarios,
''',
    '''.rc_full_gem_cache_fingerprint <- function(gem) {
  validated <- rc_validate_gem(gem)
  info <- gem$model_info %||% list()
  payload <- list(
    species = as.character(info$species %||% NA_character_),
    source = as.character(info$source %||% NA_character_),
    version = as.character(
      info$model_version %||% info$version %||% NA_character_
    ),
    commit = as.character(
      info$source_commit %||% info$commit %||% NA_character_
    ),
    checksum = as.character(info$checksum %||% NA_character_),
    S = validated$S,
    lb = validated$lb,
    ub = validated$ub
  )
  file <- tempfile("RegCompassR-gem-fingerprint-", fileext = ".rds")
  on.exit(unlink(file, force = TRUE), add = TRUE)
  saveRDS(payload, file, version = 2)
  unname(tools::md5sum(file)[[1L]])
}

#' Cache one complete full GEM per medium scenario
rc_build_full_gem_cache <- function(gem, dirs, medium_scenarios,
''',
)

replace_once(
    "R/full_gem.R",
    '''  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
''',
    '''  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  gem_fingerprint <- .rc_full_gem_cache_fingerprint(gem)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
''',
)

replace_once(
    "R/full_gem.R",
    '''    file <- file.path(
      cache_dir,
      paste0("full_gem__medium_", safe(scenario), "__condition_", safe(condition), ".rds")
    )
''',
    '''    file <- file.path(
      cache_dir,
      paste0(
        "full_gem__gem_", gem_fingerprint,
        "__medium_", safe(scenario),
        "__condition_", safe(condition), ".rds"
      )
    )
''',
)

replace_once(
    "R/full_gem.R",
    '''    if (!file.exists(file) || isTRUE(force)) {
      full <- rc_build_full_gem(
        gem = gem,
        medium_table = medium,
        condition = if (identical(condition, "all")) NULL else condition
      )
      full$condition <- condition
      saveRDS(full, file)
    } else {
      full <- readRDS(file)
    }
''',
    '''    rebuild <- !file.exists(file) || isTRUE(force)
    if (!rebuild) {
      full <- tryCatch(readRDS(file), error = function(error) NULL)
      cached_fingerprint <- if (is.list(full)) {
        full$cache_identity$gem_fingerprint %||% NA_character_
      } else {
        NA_character_
      }
      rebuild <- !identical(cached_fingerprint, gem_fingerprint)
    }
    if (rebuild) {
      full <- rc_build_full_gem(
        gem = gem,
        medium_table = medium,
        condition = if (identical(condition, "all")) NULL else condition
      )
      full$condition <- condition
      full$cache_identity <- list(
        gem_fingerprint = gem_fingerprint,
        species = gem$model_info$species %||% NA_character_,
        source = gem$model_info$source %||% NA_character_,
        version = gem$model_info$model_version %||%
          gem$model_info$version %||% NA_character_,
        commit = gem$model_info$source_commit %||%
          gem$model_info$commit %||% NA_character_,
        checksum = gem$model_info$checksum %||% NA_character_
      )
      saveRDS(full, file)
    }
''',
)

replace_once(
    "R/full_gem.R",
    '''      cache_key = paste("full_gem", scenario, condition, sep = "::"),
      medium_scenario = scenario,
''',
    '''      cache_key = paste(
        "full_gem", gem_fingerprint, scenario, condition, sep = "::"
      ),
      gem_fingerprint = gem_fingerprint,
      medium_scenario = scenario,
''',
)

# Restore positional-call compatibility while keeping species selectable by name.
replace_once(
    "R/regcompass.R",
    '''    fragment_files = NULL,
    species = c("auto", "human", "mouse"),
    sample_col = "sample_id",
''',
    '''    fragment_files = NULL,
    sample_col = "sample_id",
''',
)
replace_once(
    "R/regcompass.R",
    '''    parallel_backend = c("auto", "serial", "snow", "multicore"),
    strict_biological_defaults = TRUE,
    inference_unit = c("sample_celltype", "metacell")) {
''',
    '''    parallel_backend = c("auto", "serial", "snow", "multicore"),
    strict_biological_defaults = TRUE,
    inference_unit = c("sample_celltype", "metacell"),
    species = c("auto", "human", "mouse")) {
''',
)

replace_once(
    "R/one_shot.R",
    '''    fragment_files = NULL,
    species = c("human", "mouse"),
    gem = NULL,
    gem_version = NULL,
    humangem_version = NULL,
    medium_scenario = "physiologic",
    medium_scenarios = NULL,
    ...) {
''',
    '''    fragment_files = NULL,
    gem = NULL,
    humangem_version = NULL,
    medium_scenario = "physiologic",
    medium_scenarios = NULL,
    species = c("human", "mouse"),
    gem_version = NULL,
    ...) {
''',
)

# Forward non-default RNA assays into signed Pando confidence scoring.
replace_once(
    "R/global_workflow.R",
    '''    pando_object,
    gem,
    atac_assay = atac_assay
  )
''',
    '''    pando_object,
    gem,
    atac_assay = atac_assay,
    rna_assay = rna_assay
  )
''',
)

# Compute signed eligibility before top-k component pruning and keep module IDs sample-local.
replace_once(
    "R/pando_grn.R",
    '''    tf_peak_gene, metabolic_genes,
    top_k = top_k,
    min_shared_tfs = min_shared_tfs,
''',
    '''    tf_peak_gene, metabolic_genes,
    top_k = Inf,
    min_shared_tfs = min_shared_tfs,
''',
)

replace_once(
    "R/pando_grn.R",
    '''  relation <- as.character(edges$regulatory_relation)
  edges$used_for_component <- edges$direct_regulatory %in% TRUE |
    (!edges$direct_regulatory %in% TRUE & relation == "concordant")
  edges$used_for_component[is.na(edges$used_for_component)] <- FALSE

  for (sample in unique(as.character(nodes$sample_id))) {
''',
    '''  relation <- as.character(edges$regulatory_relation)
  component_candidate <- edges$direct_regulatory %in% TRUE |
    (!edges$direct_regulatory %in% TRUE & relation == "concordant")
  component_candidate[is.na(component_candidate)] <- FALSE
  edges$used_for_component <- component_candidate
  if (nrow(edges) && is.finite(top_k) && top_k > 0L) {
    selected <- rep(FALSE, nrow(edges))
    for (sample in unique(as.character(nodes$sample_id))) {
      sample_edge <- as.character(edges$sample_id) == sample
      sample_genes <- as.character(
        nodes$gene[as.character(nodes$sample_id) == sample]
      )
      for (gene in sample_genes) {
        index <- which(
          sample_edge & component_candidate &
            (edges$gene_a == gene | edges$gene_b == gene)
        )
        if (!length(index)) next
        order_index <- order(
          edges$direct_regulatory[index],
          edges$projection_weight[index],
          edges$shared_tf_count[index],
          decreasing = TRUE,
          na.last = TRUE
        )
        selected[index[utils::head(order_index, as.integer(top_k))]] <- TRUE
      }
    }
    edges$used_for_component <- selected
  }

  for (sample in unique(as.character(nodes$sample_id))) {
''',
)

replace_once(
    "R/pando_grn.R",
    '''    if (any(edge_index)) {
      edges$module_id[edge_index] <- nodes$module_id[
        match(
          as.character(edges$gene_a[edge_index]),
          as.character(nodes$gene)
        )
      ]
    }
''',
    '''    if (any(edge_index)) {
      sample_nodes <- nodes[node_index, , drop = FALSE]
      edges$module_id[edge_index] <- sample_nodes$module_id[
        match(
          as.character(edges$gene_a[edge_index]),
          as.character(sample_nodes$gene)
        )
      ]
    }
''',
)

test_file = r'''test_that("full-GEM cache identity changes with GEM structure", {
  make_gem <- function(stoichiometry, version) {
    S <- Matrix::Matrix(
      matrix(c(-1, stoichiometry), nrow = 1),
      sparse = TRUE,
      dimnames = list("m_e", c("EX_m", "R1"))
    )
    list(
      S = S,
      lb = stats::setNames(c(-1000, 0), colnames(S)),
      ub = stats::setNames(c(1000, 1000), colnames(S)),
      reaction_meta = data.frame(
        reaction_id = colnames(S),
        role = c("exchange", "internal"),
        stringsAsFactors = FALSE
      ),
      model_info = list(
        species = "human",
        source = "test/GEM",
        version = version,
        commit = version,
        checksum = paste0("checksum-", version)
      )
    )
  }
  medium <- data.frame(
    medium_scenario_id = "base",
    exchange_reaction_id = NA_character_,
    lb = NA_real_,
    ub = NA_real_,
    available = FALSE,
    .no_constraints = TRUE,
    stringsAsFactors = FALSE
  )
  dirs <- data.frame(
    reaction_id = "R1",
    target_direction = "forward",
    stringsAsFactors = FALSE
  )
  cache_dir <- tempfile("full-gem-cache-")
  first <- rc_build_full_gem_cache(
    make_gem(1, "v1"), dirs, medium, cache_dir = cache_dir
  )
  second <- rc_build_full_gem_cache(
    make_gem(2, "v2"), dirs, medium, cache_dir = cache_dir
  )
  first_file <- attr(first, "summary")$file[[1L]]
  second_file <- attr(second, "summary")$file[[1L]]

  expect_false(identical(first_file, second_file))
  expect_false(identical(
    attr(first, "summary")$gem_fingerprint,
    attr(second, "summary")$gem_fingerprint
  ))
  expect_equal(as.numeric(readRDS(second_file)$S["m_e", "R1"]), 2)
})

test_that("incompatible legacy species caches are invalidated", {
  S <- Matrix::Matrix(
    matrix(c(-1, 1), nrow = 1),
    sparse = TRUE,
    dimnames = list("m_e", c("EX_m", "R1"))
  )
  legacy <- list(
    S = S,
    lb = stats::setNames(c(-1000, 0), colnames(S)),
    ub = stats::setNames(c(1000, 1000), colnames(S)),
    model_info = list(source = "SysBioChalmers/Human-GEM", version = "2.0.0")
  )
  file <- tempfile(fileext = ".rds")
  saveRDS(legacy, file)
  spec <- .rc_species_gem_spec("human", "2.0.0")

  expect_warning(
    cached <- .rc_load_compatible_species_gem(file, spec),
    "Removing incompatible cached"
  )
  expect_null(cached)
  expect_false(file.exists(file))
})

test_that("species arguments preserve legacy positional ordering", {
  workflow_formals <- names(formals(rc_run_regcompass))
  one_shot_formals <- names(formals(rc_run_regcompass_one_shot))

  expect_lt(
    match("sample_col", workflow_formals),
    match("species", workflow_formals)
  )
  expect_lt(match("gem", one_shot_formals), match("species", one_shot_formals))
  expect_lt(
    match("medium_scenarios", one_shot_formals),
    match("species", one_shot_formals)
  )
})

test_that("stratum confidence forwards the selected RNA assay", {
  text <- paste(deparse(body(.rc_run_regcompass_stratum)), collapse = "\n")
  expect_match(
    text,
    "atac_assay = atac_assay,\\s*rna_assay = rna_assay"
  )
})

test_that("multi-sample projected edges receive sample-local module IDs", {
  input <- data.frame(
    sample_id = c("s1", "s1", "s2", "s2"),
    tf = c("TF1", "TF1", "TF2", "TF2"),
    target = c("A", "B", "A", "C"),
    estimate = c(1, 1, 1, 1),
    stringsAsFactors = FALSE
  )
  projected <- rc_project_metabolic_grn(
    input,
    metabolic_genes = c("A", "B", "C"),
    top_k = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    include_direct_metabolic_tf = FALSE
  )
  expect_true(all(
    startsWith(
      projected$edges$module_id,
      paste0(projected$edges$sample_id, "::")
    )
  ))
})

test_that("signed eligibility is applied before top-k component pruning", {
  input <- data.frame(
    sample_id = rep("s1", 6),
    tf = c("TF1", "TF1", "TF2", "TF2", "TF3", "TF3"),
    target = c("A", "B", "A", "C", "C", "D"),
    estimate = c(10, -10, 1, 1, 10, -10),
    stringsAsFactors = FALSE
  )
  projected <- rc_project_metabolic_grn(
    input,
    metabolic_genes = c("A", "B", "C", "D"),
    top_k = 1,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    include_direct_metabolic_tf = FALSE
  )
  pair <- paste(projected$edges$gene_a, projected$edges$gene_b, sep = "-")
  concordant <- pair == "A-C"
  discordant <- pair %in% c("A-B", "C-D")

  expect_true(any(concordant))
  expect_true(projected$edges$used_for_component[concordant])
  expect_true(all(!projected$edges$used_for_component[discordant]))
  expect_equal(
    projected$nodes$module_id[projected$nodes$gene == "A"],
    projected$nodes$module_id[projected$nodes$gene == "C"]
  )
})
'''
Path("tests/testthat/test-post-merge-audit-fixes.R").write_text(test_file, encoding="utf-8")
