#' Map reactions to metabolic modules
#' @export
rc_reaction_module_map <- function(gem, target_reactions = NULL, module_col = "metabolic_module") {
  gv <- rc_validate_gem(gem)
  meta <- gem$reaction_meta
  if (is.null(meta) || !module_col %in% colnames(meta)) {
    stop("`gem$reaction_meta` must contain module column: ", module_col, call. = FALSE)
  }
  meta <- meta[match(gv$reactions, as.character(meta$reaction_id)), , drop = FALSE]
  module <- as.character(meta[[module_col]])
  module[is.na(module) | !nzchar(module)] <- "UNASSIGNED"
  out <- data.frame(reaction_id = gv$reactions, module_id = module, stringsAsFactors = FALSE)
  if (!is.null(target_reactions)) out <- out[out$reaction_id %in% target_reactions, , drop = FALSE]
  out
}

#' Build a Module-Compass-style module meso-GEM
#' @export
rc_build_module_meso_gem <- function(gem,
                                     module_id,
                                     medium_table = NULL,
                                     condition = NULL,
                                     module_col = "metabolic_module",
                                     include_one_hop = TRUE,
                                     include_transport = TRUE,
                                     include_exchange = TRUE,
                                     include_protected = TRUE,
                                     currency_metabolites = NULL,
                                     max_reactions = 3000,
                                     strict_closure = FALSE) {
  gem <- rc_annotate_reaction_roles(gem, medium_table = medium_table)
  gv <- rc_validate_gem(gem)
  S <- gv$S
  meta <- gem$reaction_meta[match(gv$reactions, as.character(gem$reaction_meta$reaction_id)), , drop = FALSE]
  if (!module_col %in% colnames(meta)) stop("Missing module column in `reaction_meta`: ", module_col, call. = FALSE)

  module_vec <- as.character(meta[[module_col]])
  module_vec[is.na(module_vec) | !nzchar(module_vec)] <- "UNASSIGNED"
  core <- as.character(meta$reaction_id[module_vec == module_id])
  if (!length(core)) stop("No reactions found for module: ", module_id, call. = FALSE)

  keep <- core
  currency <- .rc_currency_ids(gem, currency_metabolites)
  if (isTRUE(include_one_hop)) {
    core_mets <- rownames(S)[Matrix::rowSums(abs(S[, core, drop = FALSE]) > 0) > 0]
    core_mets <- setdiff(core_mets, currency)
    if (length(core_mets)) {
      one_hop <- colnames(S)[Matrix::colSums(abs(S[core_mets, , drop = FALSE]) > 0) > 0]
      keep <- union(keep, one_hop)
    }
  }

  role <- stats::setNames(as.character(meta$role), meta$reaction_id)
  module_mets <- rownames(S)[Matrix::rowSums(abs(S[, keep, drop = FALSE]) > 0) > 0]

  touching <- function(rxns) {
    rxns <- intersect(rxns, colnames(S))
    if (!length(rxns) || !length(module_mets)) return(character())
    rxns[Matrix::colSums(abs(S[module_mets, rxns, drop = FALSE]) > 0) > 0]
  }

  if (isTRUE(include_transport)) keep <- union(keep, touching(names(role)[role == "transport"]))
  if (isTRUE(include_exchange)) {
    keep <- union(keep, touching(names(role)[role == "exchange"]))
    if (!is.null(medium_table) && "exchange_reaction_id" %in% colnames(medium_table)) {
      keep <- union(keep, as.character(medium_table$exchange_reaction_id))
    }
  }
  if (isTRUE(include_protected)) {
    protected <- names(role)[role %in% c("exchange", "transport", "demand", "sink", "maintenance", "cofactor_recycle")]
    keep <- union(keep, touching(protected))
  }

  keep <- intersect(unique(keep), colnames(S))
  max_reactions_exceeded <- is.finite(max_reactions) && length(keep) > max_reactions

  sub <- gem
  sub$S <- S[, keep, drop = FALSE]
  sub$lb <- gv$lb[keep]
  sub$ub <- gv$ub[keep]
  sub$reaction_meta <- meta[match(keep, meta$reaction_id), , drop = FALSE]
  mets_used <- rownames(sub$S)[Matrix::rowSums(abs(sub$S) > 0) > 0]
  sub$S <- sub$S[mets_used, , drop = FALSE]
  if (!is.null(gem$metabolite_meta)) {
    sub$metabolite_meta <- gem$metabolite_meta[match(mets_used, as.character(gem$metabolite_meta$metabolite_id)), , drop = FALSE]
  }
  if (!is.null(gem$gpr_table)) sub$gpr_table <- gem$gpr_table[as.character(gem$gpr_table$reaction_id) %in% keep, , drop = FALSE]

  med_diag <- data.frame()
  if (!is.null(medium_table)) {
    app <- rc_apply_medium_constraints(sub, medium_table, condition = condition, strict = FALSE)
    sub <- app$gem
    med_diag <- app$medium_diagnostics
  }
  sub$module_id <- module_id
  sub$reaction_roles <- sub$reaction_meta[, intersect(c("reaction_id", "role", "role_source", "role_confidence"), colnames(sub$reaction_meta)), drop = FALSE]
  sub$medium_diagnostics <- med_diag
  sub$closure_diagnostics <- data.frame()
  sub$build_params <- list(strategy = "module_meso_gem", module_id = module_id, max_reactions = max_reactions,
                           max_reactions_exceeded = max_reactions_exceeded,
                           n_reactions_before_any_max_guard = length(keep),
                           include_one_hop = include_one_hop, include_transport = include_transport,
                           include_exchange = include_exchange, include_protected = include_protected,
                           strict_closure = strict_closure)
  sub
}

#' Check closure diagnostics for targets in a module meso-GEM
#' @export
rc_check_module_gem_closure <- function(module_gem, target_reactions, solver = "highs", flux_threshold = 1e-8) {
  do.call(rbind, lapply(target_reactions, function(r) {
    if (!r %in% colnames(module_gem$S)) {
      return(data.frame(target_reaction = r, strict_target_feasible = FALSE, strict_vmax = NA_real_,
                        n_boundary_metabolites = NA_integer_, n_deadend_metabolites = NA_integer_,
                        n_exchange_reactions = NA_integer_, n_transport_reactions = NA_integer_,
                        n_support_reactions = NA_integer_, top_unbalanced_boundary_metabolites = NA_character_,
                        closure_warning_flag = TRUE, stringsAsFactors = FALSE))
    }
    rc_check_microgem_closure(module_gem, r, solver = solver, flux_threshold = flux_threshold)
  }))
}

#' Build module meso-GEM cache once per module and medium scenario
#' @export
rc_build_module_gem_cache <- function(gem,
                                      dirs,
                                      medium_scenarios,
                                      cache_dir = tempfile("RegCompassR_module_gem_cache_"),
                                      module_col = "metabolic_module",
                                      module_gem_params = list(),
                                      force = FALSE) {
  if (!is.data.frame(dirs) || !all(c("reaction_id", "target_direction") %in% colnames(dirs))) {
    stop("`dirs` must contain `reaction_id` and `target_direction`.", call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(medium_scenarios)) {
    medium_scenarios <- data.frame(medium_scenario_id = "base", exchange_reaction_id = character(), lb = numeric(), ub = numeric(), available = logical(), stringsAsFactors = FALSE)
  }
  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) medium_scenarios$medium_scenario_id <- "custom"

  map <- rc_reaction_module_map(gem, dirs$reaction_id, module_col = module_col)
  dirs$module_id <- map$module_id[match(dirs$reaction_id, map$reaction_id)]
  dirs$module_id[is.na(dirs$module_id) | !nzchar(dirs$module_id)] <- "UNASSIGNED"

  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  module_tasks <- unique(expand.grid(module_id = unique(dirs$module_id), medium_scenario = scenarios, stringsAsFactors = FALSE))
  module_files <- list()
  module_summary <- vector("list", nrow(module_tasks))
  for (i in seq_len(nrow(module_tasks))) {
    mid <- module_tasks$module_id[[i]]; sc <- module_tasks$medium_scenario[[i]]
    rds <- file.path(cache_dir, paste0("module_", gsub("[^A-Za-z0-9_.-]+", "_", mid), "__medium_", gsub("[^A-Za-z0-9_.-]+", "_", sc), ".rds"))
    if (!file.exists(rds) || isTRUE(force)) {
      mt <- medium_scenarios[as.character(medium_scenarios$medium_scenario_id) == sc, , drop = FALSE]
      mg <- do.call(rc_build_module_meso_gem, c(list(gem = gem, module_id = mid, medium_table = mt, module_col = module_col), module_gem_params))
      saveRDS(mg, rds)
    } else {
      mg <- readRDS(rds)
    }
    key <- paste(mid, sc, sep = "::")
    module_files[[key]] <- rds
    module_summary[[i]] <- data.frame(module_key = key, module_id = mid, medium_scenario = sc, file = rds,
                                      n_reactions = ncol(mg$S), n_metabolites = nrow(mg$S),
                                      model_version = (gem$model_info$model_version %||% gem$model_info$version %||% NA_character_),
                                      model_commit = (gem$model_info$source_commit %||% gem$model_info$commit %||% NA_character_),
                                      stringsAsFactors = FALSE)
  }

  reaction_cache <- list()
  for (i in seq_len(nrow(dirs))) {
    for (sc in scenarios) {
      row_key <- paste(dirs$reaction_id[[i]], dirs$target_direction[[i]], sc, sep = "::")
      module_key <- paste(dirs$module_id[[i]], sc, sep = "::")
      reaction_cache[[row_key]] <- list(reaction_id = dirs$reaction_id[[i]], target_direction = dirs$target_direction[[i]],
                                        medium_scenario = sc, module_id = dirs$module_id[[i]], module_key = module_key,
                                        file = module_files[[module_key]])
    }
  }
  attr(reaction_cache, "summary") <- do.call(rbind, module_summary)
  attr(reaction_cache, "target_module_map") <- dirs
  reaction_cache
}
