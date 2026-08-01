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
