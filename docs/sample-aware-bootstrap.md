# Sample-aware bootstrap contract

## Resampling rule

For condition `c`, let the observed biological samples be `D_c`. One bootstrap replicate samples `|D_c|` sample IDs with replacement from `D_c`. Every occurrence of a sampled ID contributes all cells assigned to that sample and condition. Sample IDs may be shared across conditions in paired designs because sampling is performed independently inside each condition.

This is a cluster bootstrap. It estimates reproducibility with respect to the observed donor/sample clusters rather than treating cells as independent biological replicates. Cluster sizes are preserved, so the total number of cells in a replicate can differ from the original number.

## Fallback rule

- `sample_col` supplied and valid: sample-cluster bootstrap.
- `sample_col = NULL`: print a warning and use condition-stratified cell bootstrap.
- named `sample_col` absent: print a warning naming the missing column and use condition-stratified cell bootstrap.
- existing column with missing/empty IDs: stop with an error.
- fewer than two samples in a condition: retain sample-cluster bootstrap but print a low-replication warning.

## Output provenance

Stage 1 coefficient and stability outputs include:

```text
bootstrap_method
bootstrap_resampling_unit
bootstrap_sample_col
n_bootstrap_samples_total
min_bootstrap_samples_per_condition
bootstrap_fallback_reason
```

Cell-type and condition-by-cell-type status tables additionally report sample counts. These fields are written before Stage 3 reads active edges, so downstream target-gene and core-reaction construction uses the correct stability-selected sub-GRN.

## Scope

The policy affects Stage 1 bootstrap only. Cross-validation remains condition-stratified at the cell level, Stage 2 remains condition-only, and Stage 5 still constructs one shared union GEM per medium.

The canonical multitask design does not use Pando outcome-correlation filtering. For legacy Pando sensitivity runs, the documented compatibility default remains `peak_cor = 0.01`; this value must not be inserted into the canonical `prepare_grn_design()` argument bundle.

## Timing

Each public stage prints elapsed time and final status in the R console after the final stage artifact is committed. Timing is not retained in returned objects and no timing TSV file is produced.
