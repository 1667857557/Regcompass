# Local FASTCORE and sample-balanced global calibration

## Execution phases

`rc_run_regcompass()` retains two non-overlapping process pools.

1. Each retained `condition × sample × cell-type` strict stratum is one upstream task.
2. The task constructs RNA+ATAC metacells, aggregates fragments, runs Pando, derives Pando-supported reaction confidence, expands GRN meta-modules, and completes every local meta-module with add-only FASTCORE.
3. The upstream pool is stopped and garbage-collected.
4. The all-strata/all-sample barrier verifies every retained artifact.
5. GPR-gene logCPM is combined across all metacells, optionally corrected for identifiable technical batches, and calibrated with equal total weight per biological sample.
6. Locally completed modules are unioned and deduplicated.
7. For each shared medium, the global union is checked again; FASTCORE adds support only for core directions that remain incomplete.
8. A fresh Layer 2 pool evaluates the shared model for every metacell using that metacell's own penalty vector.

Pando remains the sole peak-gene model in the integrated path. No additional Signac `LinkPeaks()` call is introduced.

## Local and global FASTCORE

Local completion is enabled by default:

```r
layer1_args = list(
  local_fastcore = TRUE,
  local_fastcore_args = list(
    solver = "highs",
    time_limit = 300,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE,
    save_models = TRUE
  )
)
```

All modules inside one strict-stratum worker reuse one unconstrained, FASTCC-screened parent GEM. This avoids rebuilding the parent for every module. Local support reactions are labeled `local_fastcore_support` and are retained when the global union is formed.

Local completion intentionally does not apply a condition-specific medium. The final global union is reconstructed once per shared, condition-invariant medium. `rc_build_meta_module_gem()` first checks every global core direction in the union and runs FASTCORE only for directions that are still incomplete.

Outputs under each stratum include:

- `03_local_fastcore/local_fastcore_models/*.rds`
- local completion summary
- directional closure diagnostics
- LP-7/LP-10 iteration diagnostics
- completed reaction membership in the stratum artifact

The merged result includes:

- `biological_reaction_membership`
- `local_completed_reaction_membership`
- `global_reaction_membership`
- `local_fastcore_summary`
- `local_fastcore_diagnostics`
- `local_fastcore_completion_iterations`

## Equal-sample global calibration

Equal-sample calibration is enabled by default:

```r
layer1_args = list(
  sample_balance = TRUE,
  sample_balance_col = "sample_id"
)
```

For sample `s` containing `n_s` metacells, each metacell receives weight:

\[
w_{su} = \frac{1/n_s}{\sum_{s'} 1}
\]

Consequently, every biological sample contributes the same total weight, regardless of the number of metacells produced. The weights are used for:

- gene-wise weighted median;
- weighted MAD and IQR robust scale;
- sigmoid gene score;
- reaction-wise weighted Q95.

GPR aggregation remains based on the fixed complete GEM GPR universe, so meta-module size does not alter gene promiscuity weights.

The global Layer 1 object records:

- `sample_balance_weights`
- `calibration_params`
- `capacity_calibration_scope`
- weighted `q95_diagnostics`

Setting `sample_balance = FALSE` restores equal-metacell calibration.

## Optional limma correction

Correction occurs after metacell logCPM merging and before gene-score calculation:

```r
layer1_args = list(
  expression_batch_correction = "limma",
  technical_batch_cols = c("library_batch"),
  preserve_design_cols = c("condition", "cell_type")
)
```

The implementation combines multiple technical columns into one interaction factor and calls `limma::removeBatchEffect()` while retaining the specified biological design.

Safety constraints:

- `sample_id` cannot be supplied as a removable technical batch;
- missing batch metadata stop the workflow;
- a single batch level produces a no-op warning;
- correction stops when batch and preserved biological design are rank-confounded;
- condition-specific medium constraints remain prohibited in shared-GEM scoring.

The Layer 1 result stores both `rna_metacell_logcpm_uncorrected` and the matrix used for capacity calculation in `rna_metacell_logcpm`, together with `expression_batch_diagnostics`.
