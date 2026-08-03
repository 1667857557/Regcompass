from pathlib import Path

path = Path('.github/scripts/refactor_celltype_union_targets.py')
text = path.read_text()
old = '''replace_once(
    stats_path,
    '      for (i in seq_len(nrow(score))) {\\n',
    '      row_indices <- which(is.na(row_meta$cell_type) |\\n'
    '        row_meta$cell_type == cell_type)\\n'
    '      for (i in row_indices) {\\n'
)
# Replace the second loop occurrence for omnibus.
text = stats_path.read_text()
old = '      for (i in seq_len(nrow(score))) {\\n'
pos = text.find(old)
if pos < 0:
    raise RuntimeError('condition statistics omnibus loop not found')
text = text[:pos] + (
    '      row_indices <- which(is.na(row_meta$cell_type) |\\n'
    '        row_meta$cell_type == cell_type)\\n'
    '      for (i in row_indices) {\\n'
) + text[pos + len(old):]
'''
new = '''text = stats_path.read_text()
old_loop = '      for (i in seq_len(nrow(score))) {\\n'
replacement_loop = (
    '      row_indices <- which(is.na(row_meta$cell_type) |\\n'
    '        row_meta$cell_type == cell_type)\\n'
    '      for (i in row_indices) {\\n'
)
if text.count(old_loop) != 2:
    raise RuntimeError(
        'condition statistics must contain exactly two scoring loops'
    )
text = text.replace(old_loop, replacement_loop, 2)
'''
if text.count(old) != 1:
    raise RuntimeError('ordered loop migration block was not found once')
path.write_text(text.replace(old, new, 1))
