from pathlib import Path

path = Path("R/condition_grn_contract.R")
text = path.read_text()
old = '''  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Install 1667857557/Pando_regcompass.", call. = FALSE)
  }
  if (!exists("infer_condition_grn", envir = asNamespace("Pando"),
              inherits = FALSE)) {
    stop("Installed Pando lacks infer_condition_grn().", call. = FALSE)
  }
'''
new = '''  pando_runtime <- .rc_require_pando_hybrid_runtime()
  .rc_step_monitor_event(
    progress_monitor, "pando_runtime",
    "validated Pando fused native condition runtime", current = 5L,
    context = list(
      pando_version = pando_runtime$version,
      native_sparse_abi = pando_runtime$native_sparse_abi,
      target_engine = pando_runtime$target_engine_backend,
      inner_cv = pando_runtime$inner_cv_backend,
      refit = pando_runtime$refit_backend
    )
  )
'''
count = text.count(old)
if count != 1:
    raise RuntimeError(f"expected one Pando runtime block, found {count}")
path.write_text(text.replace(old, new, 1))

path = Path("docs/tutorial-01-quick-start.md")
text = path.read_text()
anchor = '''```r
stopifnot(
  all(c("Group", "cell_type") %in% colnames(A@meta.data)),
  "pca" %in% names(A@reductions),
  "lsi" %in% names(A@reductions)
)
```
'''
addition = '''

Condition-aware Stage 1 requires **Pando >= 1.6.1**, native condition ABI 5,
and the registered fused C++ target engine. Inner validation uses exact
sufficient statistics with a deterministic hybrid Gram/sparse path solver;
only the outer-selected model performs per-cell OOF projection. RegCompass
checks the ABI, backend metadata and all native symbols before fitting and has
no R fallback. An incompatible installation stops immediately.
'''
count = text.count(anchor)
if count != 1:
    raise RuntimeError(f"expected one quick-start object-state anchor, found {count}")
path.write_text(text.replace(anchor, anchor + addition, 1))
