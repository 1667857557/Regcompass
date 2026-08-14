# Predefined media

`rc_make_medium_scenarios()` supports:

```text
normal_human_plasma
mouse_plasma
high_glucose
low_glucose
high_lactate
low_lactate
low_glutamine
custom
```

## Usage

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

Multiple scenarios can be requested in one call.

## Preset interpretation

- `normal_human_plasma` uses the HPLM background and retains publication provenance for encoded concentrations.
- `mouse_plasma` uses a conservative mouse-plasma availability set and mouse-specific provenance.
- glucose, lactate, and glutamine challenge presets use the same HPLM background and replace only the named challenge concentration metadata.
- human and mouse presets are species-restricted.

Published concentration values are **not** converted linearly into exchange flux bounds. They describe extracellular composition/provenance. Default exchange capacity is controlled by the medium exchange policy; an explicit `uptake_scale`, reaction-level custom bound, or other explicit flux-bound input is a modeling assumption that may change the feasible region.

Therefore two challenge presets can differ in recorded concentration while retaining the same default exchange capacity. This is intentional: concentration alone does not specify transporter kinetics or a maximum uptake flux.

## Challenge concentrations

| Scenario | Target concentration |
|---|---:|
| `high_glucose` | glucose 25 mM |
| `low_glucose` | glucose 1 mM |
| `high_lactate` | lactate 20 mM |
| `low_lactate` | lactate 0.5 mM |
| `low_glutamine` | glutamine 0.5 mM |

The output stores basal and challenge provenance separately.

## Custom reaction-level bounds

Use reaction-level `lb`/`ub` when a flux constraint is intentionally specified:

```r
custom_medium <- data.frame(
  medium_scenario_id = "my_medium",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE,
  stringsAsFactors = FALSE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

Custom bounds are intersected with the original GEM bounds; they do not open a reaction direction prohibited by the base GEM.

## Diagnostics

```r
attr(medium_scenarios, "preset_diagnostics")
attr(medium_scenarios, "medium_policy")
```

The returned table also contains scenario/reference fields that can be retained with downstream model provenance.
