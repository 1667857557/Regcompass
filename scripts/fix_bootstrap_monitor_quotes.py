from pathlib import Path
import re

path = Path("R/stepwise_workflow.R")
text = path.read_text(encoding="utf-8")
text, count = re.subn(
    r'\.rc_step_monitor_start\(\\"([^\"]+)\\", outdir, progress\)',
    r'.rc_step_monitor_start("\1", outdir, progress)',
    text,
)
if count != 6:
    raise RuntimeError(f"Expected six generated monitor lines, fixed {count}")
path.write_text(text, encoding="utf-8")
