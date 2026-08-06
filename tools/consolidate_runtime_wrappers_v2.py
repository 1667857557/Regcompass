from __future__ import annotations

import consolidate_runtime_wrappers as consolidation

_original_move_pool_wrapper = consolidation.move_pool_wrapper
_REPLACEMENT_ONLY = {
    ".rc_run_microcompass_monitored",
    ".rc_run_microcompass_engine_monitored",
}


def move_pool_wrapper(name: str, target_file: str, core_name: str) -> None:
    if name not in _REPLACEMENT_ONLY:
        _original_move_pool_wrapper(name, target_file, core_name)
        return

    # These two progress-aware functions are complete replacements: they call
    # the underlying engines directly and do not delegate to the earlier
    # minimal monitor helpers. Keep only the final implementation.
    wrapper = consolidation.function_text(consolidation.POOL, name)
    consolidation.remove_function(target_file, name)
    consolidation.remove_capture(name)
    consolidation.remove_function(consolidation.POOL, name)
    consolidation.append_wrapper(
        target_file, wrapper, f"{name} <- function"
    )


consolidation.move_pool_wrapper = move_pool_wrapper
consolidation.main()
