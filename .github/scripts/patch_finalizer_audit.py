from pathlib import Path

path = Path('.github/scripts/finalize_current_contract.py')
text = path.read_text()


def replace_once(source, old, new, label):
    count = source.count(old)
    if count != 1:
        raise RuntimeError(f'{label}: expected one match, found {count}')
    return source.replace(old, new, 1)

# Pin the current Pando PR head used by this RegCompass branch.
text = text.replace(
    "PANDO_SHA = '05246162ab5639fb12407b2b329de5149d9660a4'",
    "PANDO_SHA = '9b8cde926e176e853c824301a7daf6061d036031'",
)

# Deleted fields must not survive merely as compatibility validators.
layer1_retired = '''  retired <- c(
    "reaction_expression_condition_full_oof",
    "reaction_expression_common_oof",
    "gene_projection_condition_full_oof",
    "gene_projection_common_oof",
    "gene_projection_condition_unique_oof",
    "gene_support_condition_full_oof",
    "gene_support_common_oof",
    "gene_regulatory_modifier_condition_full_oof",
    "gene_regulatory_modifier_common_oof",
    "reaction_condition_full_support_fraction",
    "reaction_common_support_fraction"
  )
  if (any(retired %in% names(layer1))) {
    stop("Layer 1 contains retired projection routes.", call. = FALSE)
  }
'''
text = replace_once(
    text, layer1_retired, '', 'Layer 1 compatibility validator'
)
layer2_retired = '''  retired <- c(
    "penalty_condition_full_oof", "penalty_common_oof",
    "penalty_condition_unique_increment",
    "score_condition_full_oof_display_only",
    "score_common_oof_display_only",
    "score_rna_only_display_only"
  )
  if (any(retired %in% names(layer2))) {
    stop("Layer 2 contains retired penalty routes.", call. = FALSE)
  }
'''
text = replace_once(
    text, layer2_retired, '', 'Layer 2 compatibility validator'
)

# Standard-Pando fixed filters remain internal constants, not public arguments
# or output provenance fields.
standard_write = "standard_path.write_text(standard_text)\n"
standard_cleanup = '''standard_text = replace_once(
    standard_text,
    ''' + "'''" + '''      min_model_rsq = min_model_rsq,
      min_abs_estimate = max(
        .rc_standard_pando_min_abs_fixed,
        as.numeric(min_abs_estimate)
      ),
''' + "'''" + ''',
    "",
    "standard Pando obsolete threshold provenance"
)
standard_path.write_text(standard_text)
'''
text = replace_once(
    text, standard_write, standard_cleanup, 'standard Pando write hook'
)

# The current infer_condition_grn network_name is a valid naming argument.
# Only condition_grn_fit(network_name=) was retired.
text = text.replace(
    "    'network_name = \"regcompass_condition_grn\"',\n",
    "",
)

# Remove all historical route lines from current docs. Historical repository
# history remains in git; current tutorials contain only executable APIs.
text = text.replace(
    "            'penalty_common_oof', 'condition_unique', 'common_oof',\n",
    "            'condition_full_oof', 'penalty_common_oof',\n"
    "            'condition_unique', 'common_oof',\n",
)

# Replace the generated current-contract test with positive assertions only.
start_marker = "current_test.write_text(r'''"
end_marker = "''')\n\n# Static final audit."
start = text.find(start_marker)
end = text.find(end_marker, start)
if start < 0 or end < 0:
    raise RuntimeError('current contract test block was not found')
replacement = r'''current_test.write_text(r'''test_that("Layer 1 exposes one current projection route", {
  source <- paste(readLines("../../R/layer1_regulatory_support.R"),
                  collapse = "\n")
  expect_match(source, "gene_projection = projection\\$projection")
  expect_match(source, "reaction_expression = reaction_multiome")
  expect_match(source, "reaction_expression_rna_only = reaction_rna")
})

test_that("Layer 2 exposes primary and RNA-only routes", {
  source <- paste(readLines("../../R/step_layer2.R"), collapse = "\n")
  expect_match(source, 'primary = "penalty"', fixed = TRUE)
  expect_match(source, 'rna_control = "penalty_rna_only"', fixed = TRUE)
  expect_match(source, "answer$penalty_rna_only", fixed = TRUE)
})

test_that("current Pando condition API is used", {
  source <- paste(
    readLines("../../R/condition_grn_contract.R"), collapse = "\n"
  )
  expect_match(source, "Pando::condition_grn_fit(grn_object)", fixed = TRUE)
  expect_match(source, "Pando::project_condition_grn_cells", fixed = TRUE)
  expect_match(source, "Pando::infer_condition_grn", fixed = TRUE)
})

test_that("current public signatures are exact", {
  expect_identical(
    names(formals(rc_extract_pando_tf_peak_gene)),
    c("grn_object", "sample_id", "padj_threshold", "require_padj")
  )
  expect_true(all(c(
    "grn", "metacells", "meta_modules", "gem", "outdir",
    "gpr_and_method", "gene_half_saturation", "parallel", "BPPARAM",
    "progress"
  ) %in% names(formals(rc_regcompass_step_layer1))))
})
''')

# Static final audit.'''
text = text[:start] + replacement + text[end + len(end_marker):]

# Detailed diagnostics for any truly remaining retired current-contract token.
old_audit = '''found = [token for token in retired_tokens if token in source or token in doc_source]
if found:
    raise RuntimeError(f'retired current-contract tokens remain: {found}')
'''
new_audit = '''found = [token for token in retired_tokens if token in source or token in doc_source]
if found:
    locations = []
    audit_files = source_files + [value for value in text_files if value.exists()]
    for audit_path in audit_files:
        for line_number, line in enumerate(audit_path.read_text().splitlines(), 1):
            matched = [token for token in found if token in line]
            if matched:
                locations.append(
                    f"{audit_path}:{line_number}: {','.join(matched)}: {line}"
                )
    raise RuntimeError(
        'retired current-contract tokens remain:\\n' + '\\n'.join(locations)
    )
'''
if old_audit in text:
    text = replace_once(text, old_audit, new_audit, 'detailed audit')
elif 'retired current-contract tokens remain:\\n' not in text:
    raise RuntimeError('neither original nor detailed final audit was found')

path.write_text(text)
