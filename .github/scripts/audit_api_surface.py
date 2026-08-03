from __future__ import annotations

from pathlib import Path
import json
import re
import subprocess

ROOT = Path('.')
files = subprocess.check_output(['git', 'ls-files'], text=True).splitlines()
text_files = []
for name in files:
    path = ROOT / name
    if not path.is_file():
        continue
    try:
        text = path.read_text(encoding='utf-8')
    except UnicodeDecodeError:
        continue
    text_files.append((name, text))

function_pattern = re.compile(
    r'(?m)^([A-Za-z.][A-Za-z0-9._]*)\s*<-\s*function\b'
)
functions = []
for name, text in text_files:
    if not name.startswith('R/') or not name.endswith('.R'):
        continue
    for match in function_pattern.finditer(text):
        functions.append({'name': match.group(1), 'path': name,
                          'line': text.count('\n', 0, match.start()) + 1})

namespace = Path('NAMESPACE').read_text() if Path('NAMESPACE').exists() else ''
exports = re.findall(r'(?m)^export\(([^)]+)\)', namespace)
s3methods = re.findall(r'(?m)^S3method\(([^,]+),([^)]+)\)', namespace)

keywords = [
    'condition_grn', 'common_dictionary', 'common-dictionary',
    'condition_penalty', 'condition_active', 'pando', 'projection',
    'oof', 'nested_cv', 'nested-cv', 'sparse_group', 'sparse-group',
    'structural_zero', 'structural-zero', 'primary_cells', 'deprecated',
    'compatibility', 'legacy', 'alias', 'min_abs_estimate',
    'min_model_rsq', 'model_rsq', 'project_condition_grn_primary_cells',
    'project_condition_grn_to_cells', 'condition_grn_fit',
    'rc_regcompass_step_grn', 'rc_regcompass_step_layer1',
    'rc_regcompass_run', 'rc_regcompass_step_targeted_reactions'
]
reference_hits = []
for name, text in text_files:
    lower = text.lower()
    if not any(keyword.lower() in lower for keyword in keywords):
        continue
    for line_no, line in enumerate(text.splitlines(), 1):
        matched = [keyword for keyword in keywords
                   if keyword.lower() in line.lower()]
        if matched:
            reference_hits.append({
                'path': name, 'line': line_no,
                'keywords': matched, 'text': line[:500]
            })

regcompass_exports = sorted([
    value for value in exports
    if value.startswith('rc_') or any(token in value.lower() for token in
       ('condition', 'penalty', 'reaction', 'module', 'medium'))
])
relevant_functions = sorted([
    item for item in functions
    if item['name'].startswith('rc_') or any(token in item['name'].lower()
       for token in ('condition', 'dictionary', 'oof', 'structural', 'penalty'))
], key=lambda item: (item['name'], item['path'], item['line']))

report = {
    'files': files,
    'R_functions': functions,
    'relevant_functions': relevant_functions,
    'namespace_exports': exports,
    'regcompass_exports': regcompass_exports,
    's3methods': s3methods,
    'tests': [name for name in files if name.startswith('tests/')],
    'workflows': [name for name in files if name.startswith('.github/workflows/')],
    'tutorials_and_docs': [
        name for name in files
        if name.startswith(('vignettes/', 'docs/', 'man/')) or
        name.lower().startswith(('readme', 'news'))
    ],
    'reference_hits': reference_hits,
}
Path('api-surface-audit.json').write_text(
    json.dumps(report, indent=2, ensure_ascii=False)
)
with Path('api-surface-audit.txt').open('w') as handle:
    handle.write('REGCOMPASS EXPORTS\n')
    handle.write('\n'.join(regcompass_exports) + '\n\n')
    handle.write('RELEVANT FUNCTIONS\n')
    for item in relevant_functions:
        handle.write(f"{item['name']}\t{item['path']}:{item['line']}\n")
    handle.write('\nTESTS\n')
    handle.write('\n'.join(report['tests']) + '\n\n')
    handle.write('WORKFLOWS\n')
    handle.write('\n'.join(report['workflows']) + '\n\n')
    handle.write('TUTORIALS/DOCS\n')
    handle.write('\n'.join(report['tutorials_and_docs']) + '\n\n')
    handle.write('REFERENCES\n')
    for hit in reference_hits:
        handle.write(
            f"{hit['path']}:{hit['line']}\t"
            f"{','.join(hit['keywords'])}\t{hit['text']}\n"
        )
