from pathlib import Path

path = Path('tools/apply_layer2_exact_speedups.py')
text = path.read_text()
old = r'''replace_function(solver_path, '.rc_microcompass_highs_call', r'''.rc_microcompass_highs_call <- function(name, ...) {
  cache <- get0(
    ".rc_microcompass_highs_api_cache",
    envir = asNamespace("RegCompassR"), inherits = FALSE
  )
  if (!is.environment(cache)) {
    cache <- new.env(parent = emptyenv())
    assign(
      ".rc_microcompass_highs_api_cache", cache,
      envir = asNamespace("RegCompassR")
    )
  }
  fun <- get0(name, envir = cache, mode = "function", inherits = FALSE)
  if (!is.function(fun)) {
    fun <- getExportedValue("highs", name)
    assign(name, fun, envir = cache)
  }
  fun(...)
}''')'''
new = r'''replace_function(solver_path, '.rc_microcompass_highs_call', r'''.rc_microcompass_highs_call <- function(name, ...) {
  fun <- get0(
    name, envir = .rc_microcompass_highs_api_cache,
    mode = "function", inherits = FALSE
  )
  if (!is.function(fun)) {
    fun <- getExportedValue("highs", name)
    assign(name, fun, envir = .rc_microcompass_highs_api_cache)
  }
  fun(...)
}''')'''
if text.count(old) != 1:
    raise RuntimeError(f'Expected one cache block, found {text.count(old)}')
path.write_text(text.replace(old, new, 1))
print('Fixed HiGHS cache for standalone and namespace execution.')
