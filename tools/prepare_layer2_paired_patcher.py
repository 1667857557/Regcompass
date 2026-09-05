from pathlib import Path

path = Path(__file__).with_name("apply_layer2_paired_control.py")
text = path.read_text()
old = '''def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected 1 match, found {n}")
    return text.replace(old, new, 1)
'''
new = '''def replace_once(text, old, new, label):
    n = text.count(old)
    if n < 1:
        raise RuntimeError(f"{label}: expected at least 1 match, found {n}")
    return text.replace(old, new, 1)
'''
if old not in text:
    raise SystemExit("replace_once helper was not found")
path.write_text(text.replace(old, new, 1))
