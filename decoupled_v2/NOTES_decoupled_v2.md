# decoupled_v2 — updated scripts (V2 decoupled from V0)

Pushed from the local analysis session. Lives in its own directory so nothing
overlaps with existing cluster files.

## What changed vs the previous model
The repulsive wall now carries its OWN even (P2) anisotropy `zeta`, decoupled
from the dispersion tail's `a2`:

    V(R,t) = A e^{-aR} (1 + eta P1 + zeta P2) - C6 R^{-6} (1 + a1 P1 + a2 P2)

i.e. V2 is NO LONGER `a2*V0`; the lambda=2 channel now has independent radial
shape (zeta on the wall, a2 on the tail), mirroring lambda=1 (eta wall, a1 tail).
This gives the Delta j = 2 coupling genuine learnable structure.

## Generation pipeline (run on the cluster)
- `potential_generator.py` — adds `zeta` (Gaussian, sigma=0.2) as a 9th Sobol
  dimension; new `zeta` column in potentials_library.csv. Set `CONFIG["sobol_m"]`
  = 11 (2048 cand -> ~950 accepted) or 12 (4096 -> ~1900 accepted).
- `molscat_run.py` — lambda=2 wall coefficient = `zeta*A`; backward-compatible
  (`zeta = row["zeta"] if present else a2`, so old libraries still run as V2=a2*V0).
  Point `CONFIG["binary"]` at the bem MOLSCAT executable before running.

## Analysis / ML (run locally on the exported data; copies here for reference)
- `analysis/potentials.py`    — potential evaluator + DWBA overlap features.
- `analysis/logk_basis.py`    — conditioned Chebyshev basis for log k(T) in T^-1/6
                                (decorrelated regression targets; converts back to a_n).
- `analysis/features.py`, `analysis/dwba.py` — coupling features / WKB-DWBA estimator.
- `analysis/regress_kT.py`    — feature-set ablation (fingerprint vs geometry vs DWBA), CV R2.
- `analysis/analysis.py`      — fingerprint-vs-derived + dj=1 ladder + H2.
- `analysis/learning_curve.py`— R2 vs training size (data- vs information-limited).
- `analysis/hard_region.py`   — where errors and parameter-box thinness live.
- `analysis/make_valplot.py`  — out-of-fold validation figure (calibration).
- `analysis/make_examples.py`, `analysis/examples_highj.py` — example k(T) curves with GP +-1 sigma.
- `analysis/build_pdf.py`     — assembles the progress-note PDF.

Needs numpy/pandas/scipy (generation) and additionally scikit-learn (analysis).

## Headline findings so far
- Derived geometry features -> CV R2 ~ 0.90 (vs ~0.56 raw fingerprint).
- DWBA overlap adds < 0.01 (within noise); geometry is sufficient.
- GP uncertainty is well calibrated; ceiling is near-threshold physics, not data
  size (learning curve flat past ~600 potentials).
