from pathlib import Path

path = Path('.github/scripts/finalize_current_contract.py')
text = path.read_text()
text = text.replace(
    "PANDO_SHA = '05246162ab5639fb12407b2b329de5149d9660a4'",
    "PANDO_SHA = '9b8cde926e176e853c824301a7daf6061d036031'",
)
old = '''found = [token for token in retired_tokens if token in source or token in doc_source]
if found:
    raise RuntimeError(f'retired current-contract tokens remain: {found}')
'''
new = '''found = [token for token in retired_tokens if token in source or token in doc_source]
if found:
    locations = []
    audit_files = source_files + [path for path in text_files if path.exists()]
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
if text.count(old) != 1:
    raise RuntimeError(f'expected one final audit block, found {text.count(old)}')
path.write_text(text.replace(old, new, 1))
