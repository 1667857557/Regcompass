# CORDA2 in Layer 2

`rc_regcompass_step_layer2()` uses CORDA2 by default for `model_mode = "meta_module_gem"`.

Current user-facing controls are supplied through `layer2_args$model_params`, including `strict`, the CORDA2 confidence controls, and `corda2_args` (`MCxNCthresh`, `constraint`, `constrainby`, `om`, `ci`). `model_params$completion_time_limit` is not a CORDA2 control.

The current implementation treats input core reactions as an immutable structural backbone. CORDA2 runs once for each cell-type × medium structural model, and the resulting final GEM is passed directly to COMPASS-like directional scoring; there is no second post-reconstruction closure LP stage.

For equations and the exact reconstruction/scoring contract see [mathematical-model.md](mathematical-model.md). For runnable parameters see [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md).
