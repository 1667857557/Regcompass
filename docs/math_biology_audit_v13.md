# RegCompassR 1.4 mathematical and biological audit

## Scope

This audit compares the implemented workflow with FASTCORE and COMPASS and records where RegCompassR intentionally extends them. The package remains a reaction-capacity and flux-feasibility framework; it does not measure intracellular flux directly.

## Shared structural reference

The integrated workflow constructs one structural GEM per medium scenario. In `meta_module_gem` mode, the model is the union of all retained strict-stratum Pando meta-modules plus add-only FASTCORE support reactions. In `full_gem` mode, it is the complete medium-constrained Human-GEM.

All metacells share the same stoichiometric matrix, bounds, target directions and medium. Condition-specific medium rows are rejected because they would create different feasible regions.

## Evidence construction

### Expression capacity

After every retained stratum has completed, Human-GEM GPR-gene logCPM is combined across all metacells. Gene scores, GPR reaction capacities and reaction-wise Q95 calibration are recomputed on that common population. Stratum-local expression normalization is not used for cross-sample scoring.

### Pando-derived ATAC confidence

Pando is the sole peak-gene model. RegCompassR does not rerun Signac `LinkPeaks()` in the integrated workflow. Significant Pando target-region coefficients are matched to the TF-IDF peak-accessibility matrix saved by Pando. Region accessibility is transformed to a bounded metacell score and combined using absolute Pando coefficient weights. Genes without matched significant Pando regions are treated as missing ATAC evidence. Gene confidence is then aggregated through complete Human-GEM GPR groups to reaction-level confidence.

RNA capacity and Pando-derived ATAC confidence modify the objective penalty only. They do not change stoichiometry, reaction bounds or medium constraints.

## COMPASS-style directional LP

For a signed Human-GEM reaction flux vector `v`, RegCompassR first solves

\[
\max\ d v_r \quad \text{subject to}\quad Sv=0,\; l\le v\le u,
\]

where `d=1` for the forward direction and `d=-1` for reverse. It then solves

\[
\min_{v,a}\sum_j p_{j,u} a_j
\]

subject to

\[
Sv=0,\quad l\le v\le u,\quad -a_j\le v_j\le a_j,
\quad d v_r\ge \omega v_{r,\max}.
\]

This preserves the original signed Human-GEM bounds and avoids artificial forward/reverse split cycles. Penalties must match all model reactions exactly.

The evidence penalty is a RegCompass multiome extension, not the original COMPASS expression-neighbourhood formula:

\[
p_{j,u}=w_E[-\log(C_{j,u})]+w_F[-\log(F_{j,u})]
+w_M I_{missing}M+w_G m_j.
\]

Output is therefore a model-based relative metabolic potential, not an estimate of flux magnitude.

## FASTCORE completion

LP-7 identifies core reactions that can carry at least epsilon flux. LP-10 minimizes L1 flux through candidate support reactions while forcing active core directions. Reverse requests are handled by multiplying the relevant stoichiometric column by `-1` and transforming `[lb,ub]` to `[-ub,-lb]`, rather than splitting reactions.

RegCompassR uses an add-only, direction-preserving extension: all reactions selected by the global GRN/subsystem/pathway envelope are retained and FASTCORE only adds support reactions. The result is compact only with respect to reactions outside the biological envelope.

Allowed but parent-infeasible directions are reported as `parent_blocked`. Reactions with no direction allowed by their numerical bounds are reported as `no_allowed_direction`. Neither case is silently gap-filled.

## Biological and engineering constraints

- Steady-state `Sv=0` is a feasibility assumption, not proof that a reaction is active in vivo.
- Pando regulatory association does not establish causal regulation.
- Incomplete GPR AND complexes cannot define hard-core reactions or supported reaction confidence.
- Demand, sink and artificial-support reactions are not cheap default support because they can create non-biological drainage shortcuts.
- One global model is cached per medium scenario and reused for all metacells.
- The upstream strict-stratum worker pool is released before global recalibration and model construction; Layer 2 uses a fresh worker pool.
- Metacells are computational units, not biological replicates. Differential inference aggregates or models scores at the biological-sample level.

## Remaining limitations

The method has not established flux identifiability, thermodynamic direction, metabolite concentrations, enzyme kinetics or causal GRN-to-flux effects. Human-GEM annotation quality, Pando model quality and medium specification remain major determinants of feasibility. Full-scale runtime and memory benchmarks should be reported for each solver and dataset before production deployment.
