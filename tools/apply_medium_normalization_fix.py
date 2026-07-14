from pathlib import Path

path = Path("R/microcompass.R")
text = path.read_text()
old = '''  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) {
    medium_scenarios$medium_scenario_id <- "custom"
  }
  medium_scenarios$.no_constraints <- FALSE
  medium_scenarios
}'''
new = '''  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) {
    medium_scenarios$medium_scenario_id <- "custom"
  }
  if (!".no_constraints" %in% colnames(medium_scenarios)) {
    medium_scenarios$.no_constraints <- FALSE
  } else {
    medium_scenarios$.no_constraints <- as.logical(
      medium_scenarios$.no_constraints
    )
    medium_scenarios$.no_constraints[
      is.na(medium_scenarios$.no_constraints)
    ] <- FALSE
  }
  medium_scenarios
}'''
if old not in text:
    raise RuntimeError("Expected medium normalization block was not found")
path.write_text(text.replace(old, new, 1))
