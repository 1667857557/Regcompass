# Final public-contract cleanup for the global-only FASTCORE architecture.

.rc_feasibility_completion_metadata <- function(model_mode) {
  if (identical(model_mode, "meta_module_gem")) {
    return(list(
      feasibility_completion =
        "single_global_fastcore_on_each_medium_specific_union_gem",
      feasibility_completion_stages = list(
        meta_modules = paste(
          "Stage 3 defines and merges biological reaction catalogues without",
          "FASTCORE or GEM construction"
        ),
        union_gem = paste(
          "Stage 5 applies each shared medium and performs one global add-only",
          "FASTCORE completion before scoring"
        )
      )
    ))
  }
  list(
    feasibility_completion = "not_applicable_full_gem",
    feasibility_completion_stages = list(
      meta_modules =
        "Stage 3 defines complete-GPR targets but does not reconstruct a model",
      full_gem = "the complete GEM is constrained by each shared medium"
    )
  )
}

.rc_target_union_step_before_contract_cleanup <-
  rc_regcompass_step_target_union

.rc_rename_column_if_present <- function(x, old, new) {
  if (is.data.frame(x) && old %in% colnames(x)) {
    colnames(x)[colnames(x) == old] <- new
  }
  x
}

.rc_replace_exact_nonmissing <- function(x, old, new) {
  selected <- !is.na(x) & x == old
  x[selected] <- new
  x
}

.rc_clean_target_union_contract <- function(answer) {
  answer$merged_catalogue_membership <-
    answer$previous_union_membership %||% data.frame()
  answer$previous_union_membership <- NULL

  table_names <- intersect(
    c("expanded_reaction_catalog", "expanded_scoring_targets"),
    names(answer)
  )
  for (name in table_names) {
    value <- answer[[name]]
    value <- .rc_rename_column_if_present(
      value, "present_in_previous_union_membership",
      "present_in_merged_catalogue"
    )
    value <- .rc_rename_column_if_present(
      value, "previous_union_is_core", "merged_catalogue_is_core"
    )
    value <- .rc_rename_column_if_present(
      value, "previous_union_inclusion_stage",
      "merged_catalogue_inclusion_stage"
    )
    if ("target_role" %in% colnames(value)) {
      value$target_role <- .rc_replace_exact_nonmissing(
        value$target_role,
        "previous_global_core_not_rescored",
        "merged_core_not_rescored"
      )
      value$target_role <- .rc_replace_exact_nonmissing(
        value$target_role,
        "direct_database_crossref_absent_from_cached_union",
        "direct_database_crossref_absent_from_cached_union_gem"
      )
    }
    if ("lp_exclusion_reason" %in% colnames(value)) {
      value$lp_exclusion_reason <- .rc_replace_exact_nonmissing(
        value$lp_exclusion_reason,
        "absent_from_one_or_more_previous_union_models",
        "absent_from_one_or_more_cached_union_gems"
      )
    }
    answer[[name]] <- value
  }

  if (is.data.frame(answer$summary) && nrow(answer$summary)) {
    answer$summary <- .rc_rename_column_if_present(
      answer$summary,
      "n_previous_union_membership_reactions",
      "n_merged_catalogue_reactions"
    )
    if ("model_policy" %in% colnames(answer$summary)) {
      answer$summary$model_policy <-
        "reuse_exact_medium_specific_union_gems"
    }
  }

  if (is.list(answer$params)) {
    answer$params$merged_core_reactions_not_rescored <-
      answer$params$previous_core_reactions_not_rescored %||% character()
    answer$params$previous_core_reactions_not_rescored <- NULL
  }

  if (is.list(answer$microcompass)) {
    answer$microcompass$model_mode <-
      "reused_medium_specific_union_gem"
    if (is.list(answer$microcompass$params)) {
      answer$microcompass$params$shared_gem_scope <-
        "cached_medium_specific_union_gem_by_medium"
    }
    answer$microcompass$method <- paste(
      "microCOMPASS directional LP for direct",
      "KEGG/Reactome/master-Rhea-linked non-core reactions on exact cached",
      "medium-specific union GEMs"
    )
  }
  answer
}

rc_regcompass_step_target_union <- function(
    layer1, meta_modules, layer2, gem, outdir,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  answer <- .rc_target_union_step_before_contract_cleanup(
    layer1 = layer1,
    meta_modules = meta_modules,
    layer2 = layer2,
    gem = gem,
    outdir = outdir,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match,
    layer2_args = layer2_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
  answer <- .rc_clean_target_union_contract(answer)
  .rc_mm_write_tsv_gz(
    answer$merged_catalogue_membership,
    file.path(outdir, "merged_meta_module_catalogue_membership.tsv.gz")
  )
  obsolete_file <- file.path(outdir, "reused_union_gem_membership.tsv.gz")
  if (file.exists(obsolete_file)) unlink(obsolete_file)
  saveRDS(answer, file.path(outdir, "step_target_union.rds"))
  answer
}
