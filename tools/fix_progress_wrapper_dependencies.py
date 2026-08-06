from __future__ import annotations

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
STATE_REF = ".rc_layer2_progress_state$"
STATE_INIT = (
    '  progress_state <- get0(\n'
    '    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE\n'
    '  )\n'
)


def function_starts(text: str) -> list[tuple[str, int, int]]:
    pattern = re.compile(
        r"(?m)^(?P<indent>[ \t]*)(?P<name>[.A-Za-z][.A-Za-z0-9_]*)"
        r"[ \t]*<-[ \t]*function[ \t]*\("
    )
    return [
        (match.group("name"), match.start(), match.end())
        for match in pattern.finditer(text)
    ]


def find_body_open(text: str, position: int) -> int:
    quote: str | None = None
    escaped = False
    comment = False
    paren = 1
    while position < len(text):
        char = text[position]
        if comment:
            if char == "\n":
                comment = False
            position += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif char == "\\" and quote != "`":
                escaped = True
            elif char == quote:
                quote = None
            position += 1
            continue
        if char == "#":
            comment = True
        elif char in ('"', "'", "`"):
            quote = char
        elif char == "(":
            paren += 1
        elif char == ")":
            paren -= 1
            if paren == 0:
                position += 1
                break
        position += 1
    while position < len(text) and text[position].isspace():
        position += 1
    if position >= len(text) or text[position] != "{":
        raise RuntimeError("Progress wrapper must use a braced function body")
    return position


def find_body_close(text: str, body_open: int) -> int:
    quote: str | None = None
    escaped = False
    comment = False
    depth = 0
    position = body_open
    while position < len(text):
        char = text[position]
        if comment:
            if char == "\n":
                comment = False
            position += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif char == "\\" and quote != "`":
                escaped = True
            elif char == quote:
                quote = None
            position += 1
            continue
        if char == "#":
            comment = True
        elif char in ('"', "'", "`"):
            quote = char
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return position
        position += 1
    raise RuntimeError("Unclosed progress wrapper function")


def migrate_file(path: Path) -> int:
    text = path.read_text(encoding="utf-8")
    starts = function_starts(text)
    edits: list[tuple[int, int, str]] = []
    migrated = 0
    for index, (name, _start, after_formals_start) in enumerate(starts):
        next_start = starts[index + 1][1] if index + 1 < len(starts) else len(text)
        function_segment = text[after_formals_start:next_start]
        if STATE_REF not in function_segment:
            continue
        body_open = find_body_open(text, after_formals_start)
        body_close = find_body_close(text, body_open)
        body = text[body_open + 1 : body_close]
        if STATE_REF not in body:
            raise RuntimeError(
                f"State reference for {path.name}::{name} is outside its body"
            )
        if 'get0(\n    ".rc_layer2_progress_state"' in body:
            continue
        updated_body = body.replace(STATE_REF, "progress_state$")
        replacement = "{\n" + STATE_INIT + updated_body.lstrip("\n")
        edits.append((body_open, body_close, replacement))
        migrated += 1
        print(f"migrated {path.relative_to(ROOT)}::{name}")
    for start, end, replacement in reversed(edits):
        text = text[:start] + replacement + text[end:]
    if edits:
        path.write_text(text.rstrip() + "\n", encoding="utf-8")
    return migrated


def update_default() -> None:
    path = ROOT / "R/layer2_corda_pool_lifecycle.R"
    text = path.read_text(encoding="utf-8")
    old = 'model_params$model_completion %||% "fastcore"'
    new = 'model_params$model_completion %||% "corda2"'
    if old not in text:
        if new in text:
            return
        raise RuntimeError("CORDA2 request default was not found")
    path.write_text(text.replace(old, new, 1).rstrip() + "\n", encoding="utf-8")


def update_test() -> None:
    path = ROOT / "tests/testthat/test-layer2-corda-like.R"
    text = path.read_text(encoding="utf-8")
    marker = 'test_that("CORDA2 is the default structural worker-pool route"'
    if marker in text:
        return
    addition = '''

test_that("CORDA2 is the default structural worker-pool route", {
  expect_true(RegCompassR:::.rc_layer2_requested_corda2(list()))
  expect_true(RegCompassR:::.rc_layer2_requested_corda2(list(
    model_params = list()
  )))
  expect_false(RegCompassR:::.rc_layer2_requested_corda2(list(
    model_params = list(model_completion = "fastcore")
  )))
})
'''
    path.write_text(text.rstrip() + addition, encoding="utf-8")


def main() -> None:
    migrated = 0
    for path in sorted((ROOT / "R").glob("*.R")):
        if path.name == "layer2_corda_pool_lifecycle.R":
            continue
        migrated += migrate_file(path)
    if migrated == 0:
        raise RuntimeError("No progress-state-dependent wrappers were migrated")
    update_default()
    update_test()
    print(f"migrated wrappers: {migrated}")


if __name__ == "__main__":
    main()
