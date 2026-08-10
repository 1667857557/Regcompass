# Layer 2 structural model routes

`rc_regcompass_step_layer2()` supports three structural routes:

| Route | Selection | Use |
|---|---|---|
| CORDA2 | `model_mode = "meta_module_gem"` with default `model_completion` | Canonical cell-type/medium context-specific reconstruction |
| FASTCORE | `model_mode = "meta_module_gem"` and `model_params$model_completion = "fastcore"` | Supplementary context-specific reconstruction |
| Full GEM | `model_mode = "full_gem"` | Supplementary complete-network scoring |

All routes apply the selected medium by intersecting exchange bounds with the parent GEM and reuse the resulting structural model for primary and RNA-only scoring.

CORDA2 has no finite `model_params$completion_time_limit` control. FASTCORE may use completion-specific controls documented by the current code path.

See [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md) for runnable calls and [mathematical-model.md](mathematical-model.md) for the quantitative definitions.
