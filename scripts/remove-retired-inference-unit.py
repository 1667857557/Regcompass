from pathlib import Path

path = Path("R/microcompass.R")
text = path.read_text(encoding="utf-8")
old = '''  forced_unit <- getOption("RegCompassR.inference_unit", NULL)
  if (!is.null(forced_unit)) {
    unit <- match.arg(forced_unit, c("sample_celltype", "metacell"))
  }
'''
if text.count(old) != 1:
    raise RuntimeError("Expected exactly one retired inference-unit block")
path.write_text(text.replace(old, "", 1), encoding="utf-8")
