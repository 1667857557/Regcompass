from __future__ import annotations

import re

import apply_api_cleanup as cleanup


def find_function_span(
    text: str, name: str, occurrence: int = 1
) -> tuple[int, int]:
    matches = cleanup.function_matches(text, name)
    if len(matches) < occurrence:
        raise RuntimeError(
            f"Expected function {name!r} occurrence {occurrence}; found {len(matches)}"
        )
    match = matches[occurrence - 1]
    start = match.start()
    pos = match.end()

    quote: str | None = None
    escaped = False
    comment = False
    paren = 1
    while pos < len(text) and paren:
        char = text[pos]
        if comment:
            if char == "\n":
                comment = False
            pos += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif char == "\\" and quote != "`":
                escaped = True
            elif char == quote:
                quote = None
            pos += 1
            continue
        if char == "#":
            comment = True
        elif char in ('"', "'", "`"):
            quote = char
        elif char == "(":
            paren += 1
        elif char == ")":
            paren -= 1
        pos += 1
    if paren:
        raise RuntimeError(f"Unclosed formal arguments for {name!r}")

    while pos < len(text) and text[pos].isspace():
        pos += 1
    if pos >= len(text):
        raise RuntimeError(f"Missing body for {name!r}")

    if text[pos] == "{":
        depth = 0
        quote = None
        escaped = False
        comment = False
        while pos < len(text):
            char = text[pos]
            if comment:
                if char == "\n":
                    comment = False
                pos += 1
                continue
            if quote is not None:
                if escaped:
                    escaped = False
                elif char == "\\" and quote != "`":
                    escaped = True
                elif char == quote:
                    quote = None
                pos += 1
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
                    pos += 1
                    break
            pos += 1
        if depth:
            raise RuntimeError(f"Unclosed body for {name!r}")
    else:
        round_depth = square_depth = brace_depth = 0
        quote = None
        escaped = False
        comment = False
        while pos < len(text):
            char = text[pos]
            if comment:
                if char == "\n":
                    comment = False
                    if round_depth == square_depth == brace_depth == 0:
                        break
                pos += 1
                continue
            if quote is not None:
                if escaped:
                    escaped = False
                elif char == "\\" and quote != "`":
                    escaped = True
                elif char == quote:
                    quote = None
                pos += 1
                continue
            if char == "#":
                comment = True
            elif char in ('"', "'", "`"):
                quote = char
            elif char == "(":
                round_depth += 1
            elif char == ")":
                round_depth -= 1
            elif char == "[":
                square_depth += 1
            elif char == "]":
                square_depth -= 1
            elif char == "{":
                brace_depth += 1
            elif char == "}":
                brace_depth -= 1
            elif char == "\n" and round_depth == square_depth == brace_depth == 0:
                break
            pos += 1

    end = pos
    while end < len(text) and text[end] in " \t":
        end += 1
    if end < len(text) and text[end] == "\n":
        end += 1
    while end < len(text) and text[end] == "\n":
        end += 1
    return start, end


def update_description() -> None:
    path = "DESCRIPTION"
    text = cleanup.read(path)
    if "Original MATLAB CORDA2 is the default" in text:
        return
    pattern = re.compile(
        r"Stage 5 keeps\s+compact FASTCORE completion as the default and "
        r"optionally runs the original\s+MATLAB CORDA2 algorithm from "
        r"schultzdre/Constraint-Based-Modeling/CORDA2[.]m\s+after mapping "
        r"cell-type multiome evidence to HC, MC, NC and OT reactions[.]",
        re.MULTILINE,
    )
    replacement = (
        "Stage 5 uses the original MATLAB CORDA2 algorithm from "
        "schultzdre/Constraint-Based-Modeling/CORDA2.m as the default "
        "context-specific completion after mapping cell-type multiome evidence "
        "to HC, MC, NC and OT reactions. FASTCORE and complete full-GEM scoring "
        "remain explicit supplementary routes."
    )
    updated, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise RuntimeError("DESCRIPTION Layer 2 default text not found")
    # Re-wrap the DCF continuation paragraph without changing semantic content.
    lines = updated.splitlines()
    out: list[str] = []
    for line in lines:
        if not line.startswith("Description:") and not line.startswith("    "):
            out.append(line)
            continue
        if line.startswith("Description:"):
            prefix = "Description: "
            content = line[len(prefix):]
        else:
            prefix = "    "
            content = line[4:]
        if len(prefix) + len(content) <= 80:
            out.append(prefix + content)
        else:
            words = content.split()
            current = prefix
            for word in words:
                separator = "" if current.endswith(" ") else " "
                if len(current) + len(separator) + len(word) > 80:
                    out.append(current.rstrip())
                    current = "    " + word
                else:
                    current += separator + word
            out.append(current.rstrip())
    cleanup.write(path, "\n".join(out) + "\n")


cleanup.find_function_span = find_function_span
cleanup.update_description = update_description
cleanup.main()
