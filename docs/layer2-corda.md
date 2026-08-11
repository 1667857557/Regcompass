# CORDA2 Layer 2 parameters

`rc_regcompass_step_layer2()` uses CORDA2 by default when `model_mode = "meta_module_gem"`.

User-facing CORDA2 options are supplied through `layer2_args$model_params`. The current adjustable CORDA2 arguments are nested under `corda2_args`: `MCxNCthresh`, `constraint`, `constrainby`, `om`, and `ci`. `model_params$completion_time_limit` is not a CORDA2 control.

For a runnable call see [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md). For the reconstruction and scoring definitions see [mathematical-model.md](mathematical-model.md).
