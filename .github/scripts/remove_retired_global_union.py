from pathlib import Path

path = Path('R/microcompass_engine.R')
text = path.read_text()

old_guard = (
    '  mode <- match.arg(mode)\n'
    '  unit <- match.arg(unit)\n'
)
new_guard = (
    '  mode <- match.arg(mode)\n'
    '  if (!identical(mode, "full_gem")) {\n'
    '    stop("The shared full-GEM engine accepts only `full_gem` mode.",\n'
    '         call. = FALSE)\n'
    '  }\n'
    '  unit <- match.arg(unit)\n'
)
if text.count(old_guard) < 1:
    raise RuntimeError('full-GEM mode guard insertion point not found')
text = text.replace(old_guard, new_guard, 1)

start_marker = (
    '  } else {\n'
    '    model_cache <- .rc_build_medium_specific_union_gem_cache('
)
end_marker = '\n  }\n\n  units <- colnames(matrices$reaction_expression)'
start = text.find(start_marker)
end = text.find(end_marker, start)
if start < 0 or end < 0:
    raise RuntimeError('retired global union branch was not found')
replacement = (
    '  } else {\n'
    '    stop("The shared full-GEM engine received a non-full-GEM mode.",\n'
    '         call. = FALSE)\n'
    '  }'
)
text = text[:start] + replacement + text[end + len('\n  }'):]

old_module = '''          module_id = if (identical(mode, "meta_module_gem")) {
            "MEDIUM_UNION_GEM"
          } else {
            NA_character_
          },
          '''
if text.count(old_module) != 1:
    raise RuntimeError('legacy union module-ID branch was not found once')
text = text.replace(old_module, '          module_id = NA_character_,\n', 1)

old_scope = '''      shared_gem_scope = if (identical(mode, "meta_module_gem")) {
        "one_final_union_gem_per_medium_shared_across_all_units"
      } else {
        "one_full_gem_per_medium_shared_across_all_units"
      },
      '''
if text.count(old_scope) != 1:
    raise RuntimeError('legacy shared-GEM scope branch was not found once')
text = text.replace(
    old_scope,
    '      shared_gem_scope =\n'
    '        "one_full_gem_per_medium_shared_across_all_units",\n',
    1
)

method_start_marker = '    method = if (identical(mode, "full_gem")) {'
method_end_marker = '\n  )\n}\n\n#\' Run directional minimum-evidence-discordance LPs'
method_start = text.find(method_start_marker)
method_end = text.find(method_end_marker, method_start)
if method_start < 0 or method_end < 0:
    raise RuntimeError('legacy method branch structural boundary was not found')
text = (
    text[:method_start] +
    '    method = "microCOMPASS shared full-GEM directional LP"' +
    text[method_end:]
)

path.write_text(text)
