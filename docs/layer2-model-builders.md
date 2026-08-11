# Layer 2 structural routes

`rc_regcompass_step_layer2()` exposes two top-level model modes:

- `model_mode = "meta_module_gem"`: context-specific reconstruction; CORDA2 is the default completion route.
- `model_mode = "full_gem"`: complete-network scoring without context-specific reconstruction.

Within `meta_module_gem`, supplementary FASTCORE completion can be selected through the current `layer2_args$model_params` interface documented by the Rd help.

For runnable calls see [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md). Quantitative and structural definitions are maintained only in [mathematical-model.md](mathematical-model.md).
