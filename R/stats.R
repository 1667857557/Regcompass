#' Test sample-level microCOMPASS differential scores
#' @export
rc_test_microcompass_differential <- function(result, formula = score ~ condition, sample_col = "sample_id", celltype_col = "cell_type", condition_col = "condition", covariates = NULL, min_samples_per_group = 2, p_adjust_method = "BH") {
  S <- as.matrix(result$score); meta <- result$unit_meta
  rows <- lapply(rownames(S), function(r) {
    df <- data.frame(unit_id=colnames(S), score=as.numeric(S[r,]), stringsAsFactors=FALSE)
    df <- merge(df, meta, by="unit_id", all.x=TRUE)
    if (!condition_col %in% colnames(df)) return(NULL)
    ctvals <- if (celltype_col %in% colnames(df)) unique(df[[celltype_col]]) else NA_character_
    do.call(rbind, lapply(ctvals, function(ct) {
      dd <- if (!is.na(ct) && celltype_col %in% colnames(df)) df[df[[celltype_col]]==ct,,drop=FALSE] else df
      groups <- split(dd$score, dd[[condition_col]]); if (length(groups) != 2L || any(vapply(groups, length, integer(1)) < min_samples_per_group)) p <- stat <- eff <- NA_real_ else { tt <- stats::t.test(groups[[1]], groups[[2]]); p <- tt$p.value; stat <- unname(tt$statistic); eff <- mean(groups[[2]],na.rm=TRUE)-mean(groups[[1]],na.rm=TRUE) }
      data.frame(reaction_id=sub("::.*$", "", r), target_direction=sub("^.*::", "", r), cell_type=ct, contrast=paste(names(groups), collapse=" vs "), effect_size=eff, statistic=stat, p_value=p, n_samples=if (sample_col %in% colnames(dd)) length(unique(dd[[sample_col]])) else nrow(dd), n_metacells=NA_integer_, single_metacell_group_flag=NA, low_power_flag=NA, stringsAsFactors=FALSE)
    }))
  })
  out <- do.call(rbind, rows); out$FDR <- stats::p.adjust(out$p_value, method=p_adjust_method); out
}
