# Run the canonical RegCompassR 1.7.0 condition-pooled workflow.
# Cells from all samples in one condition and cell type are pooled before SuperCell2.
rc_run_regcompass <- function(
    object, gem, outdir, pfm, genome,
    fragment_files = FALSE,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    model_mode = c("meta_module_gem", "full_gem"),
    medium_scenarios = NULL,
    metacell_args = list(),
    layer1_args = list(),
    pando_args = list(),
    layer2_args = list(),
    upstream_workers = NULL,
    layer2_workers = NULL,
    parallel_backend = c("auto", "serial", "snow", "multicore"),
    strict_biological_defaults = TRUE,
    inference_unit = "metacell",
    species = c("auto", "human", "mouse")) {
  species <- .rc_infer_gem_species(gem, species)
  model_mode <- match.arg(model_mode)
  parallel_backend <- match.arg(parallel_backend)
  if (!is.character(inference_unit) || length(inference_unit) != 1L ||
      is.na(inference_unit) || !identical(inference_unit, "metacell")) {
    stop("'arg' should be one of 'metacell'", call. = FALSE)
  }
  if (!is.logical(strict_biological_defaults) ||
      length(strict_biological_defaults) != 1L ||
      is.na(strict_biological_defaults)) {
    stop("`strict_biological_defaults` must be TRUE or FALSE.", call. = FALSE)
  }
  bundles <- list(
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args,
    layer2_args = layer2_args
  )
  invalid <- names(bundles)[!vapply(bundles, is.list, logical(1))]
  if (length(invalid)) {
    stop("Workflow argument bundles must be lists: ",
         paste(invalid, collapse = ", "), call. = FALSE)
  }
  allowed_layer1 <- c(
    "regulatory_alpha", "gene_half_saturation", "tau",
    "local_fastcore", "local_fastcore_args"
  )
  unsupported_layer1 <- setdiff(names(layer1_args), allowed_layer1)
  if (length(unsupported_layer1)) {
    stop(
      "Unsupported v1.7.0 `layer1_args`: ",
      paste(unsupported_layer1, collapse = ", "),
      call. = FALSE
    )
  }
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (is.null(gem$gpr_table)) {
    stop("`gem` must contain `gpr_table`.", call. = FALSE)
  }
  validated_gem <- rc_validate_gem(gem)
  required_meta <- c(sample_col, condition_col, celltype_col)
  missing_meta <- setdiff(required_meta, colnames(object@meta.data))
  if (length(missing_meta)) {
    stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "),
         call. = FALSE)
  }
  conditions <- unique(as.character(object@meta.data[[condition_col]]))
  conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
  if (length(conditions) < 2L) {
    stop("At least two conditions are required.", call. = FALSE)
  }
  if (is.null(medium_scenarios)) {
    medium_scenarios <- rc_make_medium_scenarios(
      gem,
      scenario = "physiologic",
      species = species
    )
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(gem$model_info %||% list(), file.path(outdir, "00_model_info.rds"))
  saveRDS(medium_scenarios, file.path(outdir, "01_medium_scenarios.rds"))

  upstream_param <- .rc_phase_bpparam(upstream_workers, parallel_backend)
  on.exit(.rc_release_bpparam(upstream_param), add = TRUE)
  if (is.null(metacell_args$BPPARAM)) metacell_args$BPPARAM <- upstream_param
  pooled <- .rc_make_condition_pooled_metacells(
    object = object,
    outdir = file.path(outdir, "01_condition_pooled_metacells"),
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args,
    strict_biological_defaults = strict_biological_defaults
  )
  retained_conditions <- unique(as.character(
    pooled$metacell_meta[[condition_col]]
  ))
  missing_conditions <- setdiff(conditions, retained_conditions)
  if (length(missing_conditions)) {
    stop(
      "Condition-pooled metacell construction lost input conditions: ",
      paste(missing_conditions, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.data.frame(pooled$sample_composition) ||
      !nrow(pooled$sample_composition)) {
    stop("Condition-pooled output lacks biological-sample composition.", call. = FALSE)
  }
  metacell_object <- .rc_normalize_condition_metacell_object(
    pooled,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  if (!setequal(colnames(metacell_object), pooled$metacell_meta$metacell_id)) {
    stop("Merged metacell object and pooled metadata contain different units.",
         call. = FALSE)
  }
  saveRDS(
    metacell_object,
    file.path(outdir, "01_condition_pooled_metacells", "merged_metacell_object.rds")
  )

  meta_modules <- .rc_run_condition_pando_modules(
    metacell_object = metacell_object,
    gem = gem,
    outdir = file.path(outdir, "02_condition_pando_meta_modules"),
    pfm = pfm,
    genome = genome,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    pando_args = pando_args,
    layer1_args = layer1_args,
    BPPARAM = upstream_param
  )
  .rc_release_bpparam(upstream_param)
  upstream_param <- FALSE

  if (!is.data.frame(meta_modules$reaction_membership) ||
      !nrow(meta_modules$reaction_membership)) {
    stop("Pando meta-module construction produced no reaction membership.",
         call. = FALSE)
  }
  missing_module_reactions <- setdiff(
    unique(as.character(meta_modules$reaction_membership$reaction_id)),
    colnames(validated_gem$S)
  )
  if (length(missing_module_reactions)) {
    stop(
      "Meta-module membership contains reactions absent from the GEM: ",
      paste(utils::head(missing_module_reactions, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  global_meta_modules <- .rc_merge_stratum_meta_modules(list(list(
    group_id = "condition_pooled",
    grn_meta_modules = meta_modules
  )))
  if (!is.data.frame(global_meta_modules$global_core_reactions) ||
      !nrow(global_meta_modules$global_core_reactions)) {
    stop("No complete-GPR global core reactions remain after module merging.",
         call. = FALSE)
  }
  saveRDS(global_meta_modules, file.path(outdir, "03_global_meta_modules.rds"))

  layer1 <- .rc_build_condition_pooled_layer1(
    metacell_object = metacell_object,
    meta_modules = meta_modules,
    gem = gem,
    metacell_meta = pooled$metacell_meta,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    regulatory_alpha = layer1_args$regulatory_alpha %||% 1,
    gpr_tau = layer1_args$tau %||% 0.20,
    gene_half_saturation = layer1_args$gene_half_saturation %||%
      getOption("RegCompassR.cpm_half_saturation", 1)
  )
  if (!identical(colnames(layer1$reaction_expression),
                 as.character(layer1$unit_meta$pool_id))) {
    stop("Layer 1 reaction expression and unit metadata are not ordered identically.",
         call. = FALSE)
  }
  saveRDS(layer1, file.path(outdir, "02_global_layer1.rds"))

  cache_dir <- file.path(outdir, "04_model_cache", model_mode)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  layer2_args$model_params <- layer2_args$model_params %||% list()
  if (!is.list(layer2_args$model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  layer2_args$model_params$cache_dir <- cache_dir
  reserved_layer2 <- intersect(
    names(layer2_args),
    c(
      "layer1", "gem", "mode", "unit", "reaction_membership",
      "core_reactions", "target_reactions", "medium_scenarios",
      "sample_col", "condition_col", "celltype_col", "BPPARAM", "parallel",
      "penalty_weights"
    )
  )
  if (length(reserved_layer2)) {
    stop("`layer2_args` cannot override integrated workflow fields: ",
         paste(reserved_layer2, collapse = ", "), call. = FALSE)
  }
  layer2_param <- .rc_phase_bpparam(layer2_workers, parallel_backend)
  defaults <- list(
    layer1 = layer1,
    gem = gem,
    target_reactions = global_meta_modules$global_core_reactions$reaction_id,
    medium_scenarios = medium_scenarios,
    mode = model_mode,
    reaction_membership = if (identical(model_mode, "meta_module_gem")) {
      global_meta_modules$global_reaction_membership
    } else {
      NULL
    },
    core_reactions = if (identical(model_mode, "meta_module_gem")) {
      global_meta_modules$global_core_reactions
    } else {
      NULL
    },
    unit = "metacell",
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    parallel = TRUE,
    BPPARAM = layer2_param
  )
  defaults[names(layer2_args)] <- NULL
  microcompass <- tryCatch(
    suppressWarnings(do.call(rc_run_microcompass, c(defaults, layer2_args))),
    finally = .rc_release_bpparam(layer2_param)
  )
  if (!identical(colnames(microcompass$penalty),
                 colnames(layer1$reaction_expression))) {
    stop("microCOMPASS and Layer 1 contain different or reordered units.",
         call. = FALSE)
  }
  comparison <- .rc_condition_penalty_comparison(
    microcompass,
    condition_col = condition_col,
    celltype_col = celltype_col
  )

  result <- list(
    schema_version = "regcompass_v1.7.0_condition_pooled",
    version = "1.7.0",
    species = species,
    model_mode = model_mode,
    pooling_scope = "condition_x_celltype_across_samples",
    input_design = pooled$input_design,
    metacells = pooled,
    layer1 = layer1,
    grn_meta_modules = global_meta_modules,
    microcompass = microcompass,
    condition_summary = comparison$summary,
    condition_contrast = comparison$contrast,
    inference_policy = comparison$inference_policy,
    params = list(
      shared_gem = TRUE,
      shared_medium = TRUE,
      metacell_grouping = c(condition_col, celltype_col),
      samples_mixed_within_condition = any(
        pooled$metacell_meta$samples_mixed_within_condition
      ),
      sample_weighting = pooled$sample_weighting,
      pando_grouping = c(condition_col, celltype_col),
      inference_unit = "condition_pooled_metacell_descriptive_only",
      regulatory_alpha = layer1$capacity_params$regulatory_alpha,
      regulatory_state = "ATAC_accessibility_only",
      pando_parameter_source = "RNA_plus_ATAC_condition_x_celltype_fit",
      gpr_promiscuity_mode = "none",
      gpr_and_method = "boltzmann",
      gpr_tau = layer1$capacity_params$tau,
      gpr_or_method = "sum",
      meta_module_one_hop = TRUE,
      penalty_formula = "1/(1+log2(1+E_multiome))",
      parallel_backend = parallel_backend,
      upstream_workers = upstream_workers,
      layer2_workers = layer2_workers,
      strict_biological_defaults = strict_biological_defaults
    )
  )
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  saveRDS(result, file.path(outdir, "regcompass_condition_pooled_result.rds"))
  result
}
