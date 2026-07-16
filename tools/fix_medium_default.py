from pathlib import Path

path = Path("R/medium.R")
text = path.read_text()
old = '''    scenario = c(
      "blood_like", "minimal", "culture_like", "tumor_low_glucose",
      "low_glucose", "low_glutamine", "lactate_available", "custom"
    ),'''
new = '    scenario = "blood_like",'
if old not in text:
    raise RuntimeError("medium scenario default not found")
path.write_text(text.replace(old, new, 1))
