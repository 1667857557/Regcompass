# Give every cell-type-by-medium reconstruction an independent random stream.

.rc_complete_celltype_medium_corda_gem_noise_base <-
  .rc_complete_celltype_medium_corda_gem

.rc_corda_medium_id <- function(medium_table) {
  if (is.null(medium_table) || !is.data.frame(medium_table) ||
      !nrow(medium_table) ||
      !"medium_scenario_id" %in% colnames(medium_table)) {
    return("default")
  }
  value <- unique(as.character(medium_table$medium_scenario_id))
  value <- value[!is.na(value) & nzchar(value)]
  if (length(value) != 1L) {
    stop(
      "CORDA model construction requires exactly one medium scenario.",
      call. = FALSE
    )
  }
  value
}

.rc_corda_noise_namespace <- function(cell_type, medium_table) {
  cell_type <- as.character(cell_type)
  if (length(cell_type) != 1L || is.na(cell_type) || !nzchar(cell_type)) {
    stop("CORDA noise namespace requires one non-empty cell type.",
         call. = FALSE)
  }
  medium <- .rc_corda_medium_id(medium_table)
  paste0(
    "celltype=", utils::URLencode(cell_type, reserved = TRUE),
    "::medium=", utils::URLencode(medium, reserved = TRUE)
  )
}

.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  corda_options <- args$corda_options
  if (!is.list(corda_options)) {
    stop("CORDA options are unavailable during model construction.",
         call. = FALSE)
  }
  namespace <- .rc_corda_noise_namespace(
    cell_type = args$cell_type,
    medium_table = args$medium_table
  )
  corda_options$noise_namespace <- namespace
  args$corda_options <- corda_options
  model <- do.call(
    .rc_complete_celltype_medium_corda_gem_noise_base,
    args
  )
  model$corda_noise_contract <- list(
    seed = corda_options$seed,
    namespace = namespace,
    task_key = "namespace_x_stage_x_signed_target_x_repeat",
    distribution = paste0("Uniform(0, ", corda_options$kappa, ")"),
    scheduling_invariant = TRUE,
    matlab_rng_bitwise_identity = FALSE
  )
  model$build_params$corda_noise_namespace <- namespace
  model$build_params$corda_noise_task_key <-
    "namespace_x_stage_x_signed_target_x_repeat"
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
