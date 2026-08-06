from __future__ import annotations

from pathlib import Path

# Importing v2 executes the deterministic transformation once.
import apply_api_cleanup_v2  # noqa: F401,E402

root = Path(__file__).resolve().parents[1]
for path in list((root / "R").glob("*.R")) + [
    root / "DESCRIPTION",
    root / "README.md",
]:
    text = path.read_text(encoding="utf-8")
    path.write_text(text.rstrip() + "\n", encoding="utf-8")
