# Utility helpers will be added in later RegCompassR milestones.
NULL

.rc_first_existing_col <- function(candidates, x, fallback = NULL) {
  hit <- intersect(candidates, colnames(x))
  if (length(hit) > 0L && !is.na(hit[[1]]) && nzchar(hit[[1]])) return(hit[[1]])
  if (!is.null(fallback)) return(fallback)
  colnames(x)[1]
}
