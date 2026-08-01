from pathlib import Path

path = Path(".agent-patch/apply_detailed_progress.py")
text = path.read_text()
old = '''def insert_after_once(path, anchor, addition):
    replace_once(path, anchor, anchor + addition)
'''
new = '''def insert_after_once(path, anchor, addition):
    file = Path(path)
    current = file.read_text()
    count = current.count(anchor)
    if (
        path == "R/stepwise_workflow.R"
        and anchor == "  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)\\n"
        and count == 2
    ):
        file.write_text(current.replace(anchor, anchor + addition, 1))
        return
    replace_once(path, anchor, anchor + addition)
'''
count = text.count(old)
if count != 1:
    raise RuntimeError(
        f"expected one insert_after_once helper, found {count}"
    )
path.write_text(text.replace(old, new, 1))
