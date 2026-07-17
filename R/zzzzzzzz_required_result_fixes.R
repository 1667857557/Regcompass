# Required result-level corrections.
#
# This file is collated after the earlier compatibility layers. It changes only:
#   1. RNA library-depth normalization,
#   2. the zero/missing/no-GPR penalty ordering, and
#   3. the shared model-bound medium used by default.
#
# Pando and peak-gene links remain inferred independently within each strict
# condition x sample x cell-type stratum.

.rc_required_previous_make_supercell2_metacells <-
  rc_make_supercell2_metacells
.rc_required_previous_make_medium_scenarios <- rc_make_medium_scenarios
.rc_required_previous_run_microcompass <- rc_run_microcompass
.rc_required_previous_run_regcompass <- rc_run_regcompass
.rc_required_previous_run_regcompass_one_shot <- rc_run_regcompass_one_shot

.rc_full_library_size_cache <- new.env(parent = emptyenv())

.rc_library_cache_key <- function(ids) {
  ids <- as.character(ids)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids))) {
    stop("Metacell IDs must be non-empty before RNA normalization.",
         call. = FALSE)
  }
  paste0(length(ids), "::", paste(ids, collapse = "\001"))
}

# Cache full-transcriptome library sizes before the workflow filters RNA to GPR
# genes. Parallel workers have separate package environments, so each stratum
# stores and consumes its own entry.
rc_make_supercell2_metacells <- function(...) {
  answer <- .rc_required_previous_make_supercell2_metacells(...)
  counts <- answer$rna_counts
  if (!is.null(counts) && !is.null(dim(counts)) && ncol(counts) > 0L) {
    ids <- colnames(counts)
    key <- .rc_library_cache_key(ids)
    library_size <- Matrix::colSums(counts)
    if (any(!is.finite(library_size)) || any(library_size <= 0)) {
      stop("Every metacell must have a positive finite full RNA library size.",
           call. = FALSE)
    }
    assign(key, library_size, envir = .rc_full_library_size_cache)
  }
  answer
}

.rc_metacell_logcpm <- function(counts, scale_factor = 1e6,
                                library_size = NULL) {
  counts <- methods::as(counts, "dgCMatrix")
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("`scale_factor` must be one positive finite number.", call. = FALSE)
  }

  normalization_scope <- "input_matrix_library_size"
  if (is.null(library_size)) {
    key <- .rc_library_cache_key(colnames(counts))
    if (exists(key, envir = .rc_full_library_size_cache, inherits = FALSE)) {
      library_size <- get(
        key,
        envir = .rc_full_library_size_cache,
        inherits = FALSE
      )
      rm(list = key, envir = .rc_full_library_size_cache)
      normalization_scope <-
        "full_transcriptome_library_size_before_gpr_filter"
    } else {
      library_size <- Matrix::colSums(counts)
    }
  }

  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts)) {
    stop("`library_size` must contain one value per metacell.",
         call. = FALSE)
  }
  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop("Every metacell must have a positive finite RNA library size.",
         call. = FALSE)
  }

  scaled <- counts %*% Matrix::Diagonal(
    x = scale_factor / library_size
  )
  answer <- log1p(scaled)
  attr(answer, "normalization_scope") <- normalization_scope
  attr(answer, "library_size") <- stats::setNames(
    library_size,
    colnames(counts)
  )
  answer
}

.rc_compass_model_bound_medium <- function(gem, exchange_limit = 1) {
  if (!is.numeric(exchange_limit) || length(exchange_limit) != 1L ||
      !is.finite(exchange_limit) || exchange_limit <= 0) {
    stop("`exchange_limit` must be one positive finite number.",
         call. = FALSE)
  }

  validated <- rc_validate_gem(gem)
  if (is.null(gem$reaction_meta) ||
      !"role" %in% colnames(gem$reaction_meta)) {
    gem <- rc_annotate_reaction_roles(gem)
  }
  meta <- gem$reaction_meta[
    match(validated$reactions, as.character(gem$reaction_meta$reaction_id)),
    ,
    drop = FALSE
  ]
  exchange <- as.character(meta$reaction_id[
    as.character(meta$role) == "exchange"
  ])
  exchange <- intersect(exchange, validated$reactions)
  if (!length(exchange)) {
    stop(
      "No exchange reactions were identified for COMPASS-style model bounds.",
      call. = FALSE
    )
  }

  index <- match(exchange, validated$reactions)
  original_lb <- as.numeric(validated$lb[index])
  original_ub <- as.numeric(validated$ub[index])
  lb <- pmax(original_lb, -exchange_limit)
  ub <- pmin(original_ub, exchange_limit)
  if (any(lb > ub)) {
    stop("Capping exchange bounds produced invalid lower/upper bounds.",
         call. = FALSE)
  }

  data.frame(
    medium_scenario_id = "compass_model_bounds",
    exchange_reaction_id = exchange,
    metabolite_id = if ("metabolite_id" %in% colnames(meta)) {
      as.character(meta$metabolite_id[index])
    } else {
      NA_character_
    },
    condition = "all",
    lb = lb,
    ub = ub,
    available = TRUE,
    original_lb = original_lb,
    original_ub = original_ub,
    exchange_limit = exchange_limit,
    evidence_source =
      "gem_directionality_with_compass_uniform_exchange_limit",
    assumption_level = "shared_model_defined_environment",
    target_exchange_flag = FALSE,
    concentration_used_for_rate_bound = FALSE,
    rate_bound_source = "gem_bounds_capped_like_compass",
    stringsAsFactors = FALSE
  )
}

rc_make_medium_scenarios <- function(
    gem,
    scenario = "compass_model_bounds",
    custom_medium = NULL,
    uptake_scale = c(
      permissive_all_exchange = 1,
      blood_like = 1, culture_like = 1, minimal = 0.1,
      tumor_low_glucose = 0.5, low_glucose = 0.1,
      low_glutamine = 0.1, lactate_available = 1
    ),
    condition_col = NULL,
    exchange_roles = c("exchange"),
    condition = condition_col,
    exchange_limit = 1) {
  choices <- c(
    "compass_model_bounds",
    "permissive_all_exchange", "blood_like", "culture_like", "minimal",
    "tumor_low_glucose", "low_glucose", "low_glutamine",
    "lactate_available", "custom"
  )
  scenario <- match.arg(scenario, choices = choices, several.ok = TRUE)
  pieces <- list()

  if ("compass_model_bounds" %in% scenario) {
    condition_value <- as.character(condition %||% "all")
    if (length(condition_value) != 1L || is.na(condition_value) ||
        !condition_value %in% c("", "all")) {
      stop(
        "COMPASS-style model bounds must be shared across conditions.",
        call. = FALSE
      )
    }
    pieces[[length(pieces) + 1L]] <- .rc_compass_model_bound_medium(
      gem,
      exchange_limit = exchange_limit
    )
  }

  other <- setdiff(scenario, "compass_model_bounds")
  if (length(other)) {
    pieces[[length(pieces) + 1L]] <-
      .rc_required_previous_make_medium_scenarios(
        gem = gem,
        scenario = other,
        custom_medium = custom_medium,
        uptake_scale = uptake_scale,
        condition_col = condition_col,
        exchange_roles = exchange_roles,
        condition = condition
      )
  }

  .rc_bind_frames_fill(pieces)
}

# COMPASS maps missing reaction expression to zero before applying an inverse
# expression penalty. Since C_abs = x / (x + 1), 1 - C_abs = 1 / (x + 1).
# Observed zero, missing expression, and no-GPR reactions therefore receive the
# same maximum expression penalty instead of rewarding missing annotation.
.rc_compute_multiome_penalty_core <- function(
    C_rel, reaction_confidence, gpr_diagnostics = NULL,
    reaction_roles = NULL,
    weights = c(expr = 1.0, confidence = 0.5, missing = 1.0,
                gpr_missing = 0),
    eps = 1e-6, penalty_cap = 20,
    support_penalty = c(
      exchange = 1.0, demand = 20, sink = 20,
      artificial_support = 20, cofactor_recycle = 0.50, transport = 1.00
    ),
    missing_penalty = 1) {
  if (!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps <= 0 ||
      !is.numeric(penalty_cap) || length(penalty_cap) != 1L ||
      !is.finite(penalty_cap) || penalty_cap <= 0 ||
      !is.numeric(missing_penalty) || length(missing_penalty) != 1L ||
      !is.finite(missing_penalty) || missing_penalty < 0) {
    stop(
      "Penalty constants must be finite and satisfy eps > 0, cap > 0, missing >= 0.",
      call. = FALSE
    )
  }

  C_input <- as.matrix(C_rel)
  F_input <- rc_layer2_confidence_matrix(
    reaction_confidence,
    C_input
  )
  if (is.null(rownames(C_input)) || is.null(colnames(C_input)) ||
      is.null(rownames(F_input)) || is.null(colnames(F_input)) ||
      anyDuplicated(rownames(C_input)) ||
      anyDuplicated(colnames(C_input)) ||
      anyDuplicated(rownames(F_input)) ||
      anyDuplicated(colnames(F_input))) {
    stop(
      "Capacity and confidence matrices require unique reaction and unit IDs.",
      call. = FALSE
    )
  }

  reactions <- union(rownames(C_input), rownames(F_input))
  units <- union(colnames(C_input), colnames(F_input))
  align <- function(matrix_in, fill = NA_real_) {
    output <- matrix(
      fill,
      nrow = length(reactions),
      ncol = length(units),
      dimnames = list(reactions, units)
    )
    common_r <- intersect(reactions, rownames(matrix_in))
    common_u <- intersect(units, colnames(matrix_in))
    output[common_r, common_u] <-
      matrix_in[common_r, common_u, drop = FALSE]
    output
  }

  C <- align(C_input)
  F_original <- align(F_input)
  C[is.finite(C)] <- .rc_clamp01(C[is.finite(C)])
  F_original[is.finite(F_original)] <-
    .rc_clamp01(F_original[is.finite(F_original)])

  missing_expression_flag <- !is.finite(C)
  P_expr <- 1 - C
  P_expr[missing_expression_flag] <- missing_penalty

  finite_regulation <- is.finite(F_original)
  F_effective <- matrix(
    1,
    nrow = nrow(F_original),
    ncol = ncol(F_original),
    dimnames = dimnames(F_original)
  )
  F_effective[finite_regulation] <- pmin(
    1,
    pmax(2 * F_original[finite_regulation], eps)
  )
  P_conf <- -log(pmax(F_effective, eps))
  P_conf[!is.finite(P_conf)] <- 0

  # Missing expression already receives the same maximum inverse-expression
  # penalty as observed zero; do not add a second missingness penalty.
  P_missing <- matrix(
    0,
    nrow = nrow(C),
    ncol = ncol(C),
    dimnames = dimnames(C)
  )

  P_gpr <- matrix(
    0,
    nrow = nrow(C),
    ncol = ncol(C),
    dimnames = dimnames(C)
  )
  gpr_missing_fraction <- stats::setNames(
    rep(0, nrow(C)),
    rownames(C)
  )
  if (!is.null(gpr_diagnostics)) {
    if (!is.data.frame(gpr_diagnostics) ||
        !all(c("reaction_id", "missing_gene_fraction") %in%
             colnames(gpr_diagnostics))) {
      stop(
        "`gpr_diagnostics` must contain reaction_id and missing_gene_fraction.",
        call. = FALSE
      )
    }
    if (anyDuplicated(as.character(gpr_diagnostics$reaction_id))) {
      stop("`gpr_diagnostics$reaction_id` must be unique.",
           call. = FALSE)
    }
    hit <- intersect(
      rownames(C),
      as.character(gpr_diagnostics$reaction_id)
    )
    values <- as.numeric(
      gpr_diagnostics$missing_gene_fraction[
        match(hit, as.character(gpr_diagnostics$reaction_id))
      ]
    )
    values[!is.finite(values)] <- 0
    gpr_missing_fraction[hit] <- pmin(pmax(values, 0), 1)
    P_gpr <- matrix(
      gpr_missing_fraction[rownames(C)],
      nrow = nrow(C),
      ncol = ncol(C),
      dimnames = dimnames(C)
    )
  }

  default_weights <- c(
    expr = 1,
    confidence = 0.5,
    missing = 1,
    gpr_missing = 0
  )
  if (is.null(names(weights)) ||
      any(!names(weights) %in% names(default_weights)) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop(
      "`weights` must be a named non-negative vector using expr, confidence, missing, or gpr_missing.",
      call. = FALSE
    )
  }
  W <- default_weights
  W[names(weights)] <- weights
  P_base <- W[["expr"]] * P_expr +
    W[["confidence"]] * P_conf +
    W[["missing"]] * P_missing +
    W[["gpr_missing"]] * P_gpr

  role <- stats::setNames(rep("internal", nrow(C)), rownames(C))
  role_source <- stats::setNames(
    rep("unknown", nrow(C)),
    rownames(C)
  )
  role_confidence <- stats::setNames(
    rep(NA_character_, nrow(C)),
    rownames(C)
  )
  if (!is.null(reaction_roles)) {
    roles <- if (is.data.frame(reaction_roles)) {
      reaction_roles
    } else {
      as.data.frame(reaction_roles)
    }
    if (!all(c("reaction_id", "role") %in% colnames(roles))) {
      stop(
        "`reaction_roles` must contain reaction_id and role.",
        call. = FALSE
      )
    }
    if (anyDuplicated(as.character(roles$reaction_id))) {
      stop("`reaction_roles$reaction_id` must be unique.",
           call. = FALSE)
    }
    hit <- intersect(rownames(C), as.character(roles$reaction_id))
    index <- match(hit, as.character(roles$reaction_id))
    role[hit] <- as.character(roles$role[index])
    if ("role_source" %in% colnames(roles)) {
      role_source[hit] <- as.character(roles$role_source[index])
    }
    if ("role_confidence" %in% colnames(roles)) {
      role_confidence[hit] <- as.character(roles$role_confidence[index])
    }
  }

  if (is.null(names(support_penalty)) ||
      any(!is.finite(support_penalty)) ||
      any(support_penalty < 0)) {
    stop(
      "`support_penalty` must be a named finite non-negative vector.",
      call. = FALSE
    )
  }

  structural_role <- role %in%
    c("exchange", "demand", "sink", "artificial_support")
  curated_role <- role %in% names(support_penalty) &
    role_source %in% c("curated", "model_high_confidence")
  override <- (structural_role | curated_role) &
    role %in% names(support_penalty)
  support_penalty_used <- stats::setNames(
    rep(NA_real_, nrow(C)),
    rownames(C)
  )
  support_penalty_used[override] <- as.numeric(
    support_penalty[role[override]]
  )

  P <- P_base
  if (any(override)) {
    P[override, ] <- support_penalty_used[override]
  }
  P <- pmin(pmax(P, 0), penalty_cap)
  P[!is.finite(P)] <- penalty_cap

  P_role <- matrix(
    0,
    nrow = nrow(C),
    ncol = ncol(C),
    dimnames = dimnames(C)
  )
  P_role[override, ] <- P[override, , drop = FALSE]

  list(
    penalty = P,
    components = list(
      P_expr = P_expr,
      P_conf = P_conf,
      P_missing = P_missing,
      P_gpr = P_gpr,
      P_role = P_role,
      P_base = P_base,
      C_abs = C,
      C_rel = C,
      reaction_regulatory_support = F_original,
      reaction_confidence = F_original,
      reaction_confidence_effective = F_effective,
      missing_expression_flag = missing_expression_flag,
      missing_regulatory_support_flag = !finite_regulation,
      missing_evidence_flag =
        missing_expression_flag | !finite_regulation,
      gpr_missing_fraction = gpr_missing_fraction,
      role = role,
      role_source = role_source,
      role_confidence = role_confidence,
      role_override_flag = override,
      support_penalty_used = support_penalty_used
    ),
    evidence_policy = paste(
      "COMPASS-like inverse-expression penalty from bounded absolute RNA",
      "support; observed zero, missing expression and no-GPR reactions share",
      "the maximum expression penalty; exchange flux is controlled by one",
      "shared model-bound medium"
    ),
    penalty_formula = paste(
      "w_expr*(1-C_abs) + repression_modifier +",
      "w_gpr*gpr_missing; missing C_abs uses the same maximum",
      "inverse-expression penalty as C_abs=0"
    )
  )
}

rc_run_microcompass <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("sample_celltype", "metacell"),
    condition_col = "condition", sample_col = "sample_id",
    celltype_col = "cell_type", model_params = list(),
    penalty_weights = c(expr = 1.0, confidence = 0.5, missing = 1.0),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    time_limit = 60, flux_threshold = 1e-8,
    BPPARAM = NULL) {
  if (is.null(medium_scenarios) && is.null(medium_table)) {
    medium_scenarios <- rc_make_medium_scenarios(
      gem,
      scenario = "compass_model_bounds"
    )
  }
  .rc_required_previous_run_microcompass(
    layer1 = layer1,
    gem = gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = mode,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    unit = unit,
    condition_col = condition_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    model_params = model_params,
    penalty_weights = penalty_weights,
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM
  )
}

rc_run_regcompass <- function(
    object, gem, outdir, pfm, genome,
    fragment_files = NULL,
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
    inference_unit = c("sample_celltype", "metacell")) {
  if (is.null(medium_scenarios)) {
    medium_scenarios <- rc_make_medium_scenarios(
      gem,
      scenario = "compass_model_bounds"
    )
  }
  answer <- .rc_required_previous_run_regcompass(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    model_mode = model_mode,
    medium_scenarios = medium_scenarios,
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args,
    layer2_args = layer2_args,
    upstream_workers = upstream_workers,
    layer2_workers = layer2_workers,
    parallel_backend = parallel_backend,
    strict_biological_defaults = strict_biological_defaults,
    inference_unit = inference_unit
  )
  answer$params$medium_policy <-
    "shared_compass_model_bounds_unless_explicitly_overridden"
  answer
}

rc_run_regcompass_one_shot <- function(
    object, outdir, pfm, genome,
    fragment_files = NULL,
    gem = NULL,
    humangem_version = "2.0.0",
    medium_scenario = "compass_model_bounds",
    medium_scenarios = NULL,
    ...) {
  .rc_required_previous_run_regcompass_one_shot(
    object = object,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    gem = gem,
    humangem_version = humangem_version,
    medium_scenario = medium_scenario,
    medium_scenarios = medium_scenarios,
    ...
  )
}
