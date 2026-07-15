# RegCompassR v1.3 mathematical and biological audit

## Scope

This audit compares the implemented workflow with the original FASTCORE and
COMPASS algorithms and records where RegCompassR intentionally extends them.
The package remains a reaction-capacity and flux-feasibility framework; it does
not measure intracellular flux directly.

## COMPASS-style directional LP

For a signed Human-GEM reaction flux vector `v`, RegCompassR first solves

\[
\max\ d v_r \quad \text{subject to}\quad Sv=0,\; l\le v\le u,
\]

where `d=1` for the forward direction and `d=-1` for reverse. It then solves

\[
\min_{v,a}\sum_j p_j a_j
\]

subject to

\[
Sv=0,\quad l\le v\le u,\quad -a_j\le v_j\le a_j,
\quad d v_r\ge \omega v_{r,\max}.
\]

This is equivalent to penalizing two non-negative directional copies when both
copies use the same penalty, but it preserves the original signed Human-GEM
bounds and avoids artificial forward/reverse split cycles. Penalties must match
all model reactions exactly; missing penalties are now errors rather than being
silently replaced.

The evidence penalty is a RegCompass multiome extension, not the original
COMPASS expression/neighbourhood formula:

\[
p_{j,u}=w_E[-\log(C_{j,u})]+w_F[-\log(F_{j,u})]
+w_M I_{missing}M+w_G m_j.
\]

Accordingly, output should be interpreted as a model-based relative metabolic
potential, not an estimate of flux magnitude.

## FASTCORE completion

The LP-7 and LP-10 implementations retain the original roles: LP-7 identifies
core reactions that can carry at least epsilon flux, and LP-10 minimizes the
L1 flux carried by non-core reactions while forcing the active core. Reverse
requests are handled by multiplying the relevant stoichiometric column by -1
and transforming `[lb,ub]` to `[-ub,-lb]`, rather than splitting reactions.

RegCompassR intentionally uses an **add-only, direction-preserving extension**:
all GRN/subsystem/pathway biological members are retained and FASTCORE only
adds support reactions. With `target_direction="both"`, every parent-feasible
direction is preserved. This is stricter and generally larger than canonical
FASTCORE, which only requires each core reaction to be active in at least one
feasible mode. The result is therefore compact only with respect to reactions
outside the biological envelope.

## Biological and engineering constraints

- Steady-state `Sv=0` is a feasibility assumption, not proof that a reaction is
  active in vivo.
- RNA, ATAC, GRN and GPR evidence alter penalties or biological membership; they
  do not override stoichiometry or medium bounds.
- Incomplete GPR AND complexes cannot define hard core reactions.
- Demand, sink and artificial-support reactions are not cheap default support;
  otherwise they can create non-biological drainage shortcuts.
- Medium-constrained parent models are cached by scenario and condition.
- Metacells are computational units, not biological replicates. Differential
  inference aggregates scores to one value per biological sample and cell type.
- A blocked hard core is reported as `no_allowed_direction`; it is never silently
  reported as a successful reconstruction.

## Remaining limitations

The method has not established flux identifiability, thermodynamic direction,
metabolite concentrations, enzyme kinetics, or causal GRN-to-flux effects.
Human-GEM annotation quality and medium specification remain major determinants
of feasibility. Full-scale Human-GEM runtime and memory benchmarks should be
reported for each solver and dataset before production deployment.
