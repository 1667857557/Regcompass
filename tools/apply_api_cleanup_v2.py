from __future__ import annotations

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


cleanup.find_function_span = find_function_span
cleanup.main()
