# Fictitious Potential Generator — Implementation Plan

Standalone generator (independent of MOLSCAT) that Sobol-samples a family of
simple atom + linear-rigid-rotor potentials, validates each one physically, and
writes a CSV library of MOLSCAT-ready coefficients plus characterization
descriptors.

---

## 1. Scope and goal

- **Target collision type:** MOLSCAT `ITYPE = 1` (atom + linear rigid rotor).
- **Potential model:** isotropic V0 plus two leading anisotropies (V1, V2),
  expanded in **unnormalized** Legendre polynomials (MOLSCAT convention:
  `P1 = cos θ`, `P2 = 3 cos²θ − 1`).
- **Isotropic term (Born–Mayer + dispersion, no damping):**

  ```
  V0(R) = A · exp(−a · R) − C6 · R^(−6)
  ```

- **Anisotropic terms as scalings of V0:**

  ```
  V1(R) = a1 · V0(R)
  V2(R) = a2 · V0(R)
  ```

- **Full angular potential:**

  ```
  V(R, θ) = V0(R) · [ 1 + a1·P1(cosθ) + a2·P2(cosθ) ]
  ```

- **Output:** a validated CSV library — one row per surviving potential — with
  raw parameters, derived well descriptors at θ = 0°, 90°, 180°, anisotropy
  contrasts, and the per-system rotor/mass data (B, M, μ).

**Design note:** no damping function is used (to keep MOLSCAT input simple), so
the `−C6 R^(−6)` term diverges at small R. The generator must therefore
**screen out inner-collapse** candidates rather than rely on a damping term.

---

## 2. Parameter space and sampling

| Parameter | Symbol | Scale | Role |
|-----------|--------|-------|------|
| Repulsive prefactor | A | log-uniform | V0 repulsion strength |
| Repulsive exponent | a | log-uniform | V0 repulsive steepness |
| Dispersion coefficient | C6 | log-uniform | V0 long-range attraction |
| V1 scaling | a1 | linear, symmetric about 0 | odd anisotropy (heteronuclear only) |
| V2 scaling | a2 | linear, symmetric about 0 | quadrupolar anisotropy |
| Rotational constant | B | log-uniform | rotor level spacing |
| Molecular mass | M | log-uniform | → reduced mass μ |

### Sampling method
- Use a **Sobol sequence** (`scipy.stats.qmc.Sobol`) over the 7-D unit cube,
  then scale each dimension to its physical range.
- Sobol gives better space-filling than uniform random at every sample count and
  extends cleanly when later appending active-learning batches.
- Request `2^m` points for Sobol balance properties (e.g. 256, 1024, 4096).

### Range definition
- Ranges are **config constants** at the top of the script (easy to retune).
- `A`, `a`, `C6` ranges chosen so the *typical* well depth/position lands in a
  realistic van der Waals regime — but final coverage is checked in **derived**
  (ε, R_m) space, not raw space (see §6).
- `|a1|, |a2|` bounded (e.g. ≤ 0.3–0.5) so the angular factor stays physical.
- `B ∈ [1e-3, 10] cm⁻¹`; `M` over the target mass range.

### Reduced mass
- Sample **M**, then compute μ = m·M / (m + M) with m = 2 (H₂) or m = 4 (He).
- Never sample μ directly (it saturates near m for heavy M).

### Optional B–M coupling
- Heavy molecules physically have small B. Independent sampling produces
  unphysical heavy-but-stiff rotors.
- Optional: impose a rough trend `B̄ ∝ M^(−α)` with scatter, or reject
  implausible (B, M) corners. Off by default; flag in config.

---

## 3. Per-candidate construction

For a sampled `(A, a, C6, a1, a2)`:

1. Define `V0(R) = A·exp(−a·R) − C6·R^(−6)`.
2. Define the angular factor `g(θ) = 1 + a1·cosθ + a2·(3cos²θ − 1)`.
3. The radial cut at fixed θ is `V(R, θ) = V0(R) · g(θ)` — i.e. a scaled V0.

Because each angular cut is just a scaled V0, characterization at the three
special angles reduces to analyzing V0 scaled by `g(0)`, `g(90°)`, `g(180°)`.

Angular factor values (unnormalized Legendre):
- θ = 0°   : `g = 1 + a1·(1) + a2·(2)`   ( cos0 = 1,  3·1−1 = 2 )
- θ = 90°  : `g = 1 + a1·(0) + a2·(−1)`  ( cos90 = 0, 3·0−1 = −1 )
- θ = 180° : `g = 1 + a1·(−1) + a2·(2)`  ( cos180 = −1, 3·1−1 = 2 )

---

## 4. Characterization (per surviving potential)

For each angle θ ∈ {0°, 90°, 180°}:

- **Outer minimum:** solve `dV/dR = 0` for the outer root → well depth `ε(θ)`
  and minimum position `R_m(θ)`.
  - `dV0/dR = −a·A·exp(−a·R) + 6·C6·R^(−7)`.
  - Find sign change of dV/dR on a radial grid, then refine (Brent's method).
- **Repulsive-wall radius:** the inner zero-crossing `V(R, θ) = 0`
  (energy-independent descriptor; the σ-like radius).
  *(Default choice — see §7. Inner classical turning point at a reference
  energy is an alternative if requested.)*

### Derived anisotropy descriptors
- `Δε_odd  = ε(0°) − ε(180°)`  → odd (V1-like) asymmetry.
- `Δε_even = ½·(ε(0°) + ε(180°)) − ε(90°)`  → quadrupolar (V2-like) anisotropy.
- These are physically interpretable and serve as good ML features and
  stage-2 mapping targets.

---

## 5. Screening / rejection logic

A candidate is **accepted** only if all hold:

1. **Single outer well** at each of θ = 0°, 90°, 180° (exactly one outer
   minimum on the radial grid).
2. **No inner collapse:** scanning inward from `R_m(θ)`, V rises monotonically
   to `R_min` (the `−C6 R^(−6)` term does not overtake the exponential before
   the wall). Reject if V turns back down.
3. **Angular physicality (full scan):** `g(θ) > g_floor` for all
   θ ∈ [0, π] (cheap 1-D scan in cosθ). Guards against the well inverting or
   vanishing at *intermediate* angles even when 0/90/180 look fine — the one
   thing three-angle sampling can miss. Use a small positive floor, not just 0.
4. **Physical well depth:** `ε(θ)` within sane bounds (not absurdly deep/shallow).

Record a **rejection reason** for every discarded candidate (provenance) to
understand coverage and tune ranges.

---

## 6. Coverage check (post-generation)

- Uniform Sobol coverage in `(A, a, C6)` does **not** imply uniform coverage in
  the rate-relevant `(ε, R_m)` space.
- After generation: histogram / 2-D scatter of derived `(ε(90°), R_m(90°))`,
  `Δε_odd`, `Δε_even`. Inspect for gaps or clustering.
- If coverage is poor, either widen/retune raw ranges or reweight/resample
  toward under-covered (ε, R_m) regions.

---

## 7. Decisions to confirm (defaults used unless changed)

1. **Inner turning point** = `V(R) = 0` crossing (energy-independent).
   Alternative: classical turning point at a specified reference energy.
2. **Legendre convention** = **unnormalized** (`P2 = 3cos²θ − 1`), matching
   MOLSCAT, so the generated/screened potential equals what MOLSCAT integrates.
3. **B–M coupling** = independent by default; optional trend available.
4. **Damping** = none (by design).

---

## 8. CSV output schema

One row per **accepted** potential:

| Column | Description |
|--------|-------------|
| `id` | unique index |
| `A`, `a`, `C6` | raw V0 (Born–Mayer + dispersion) parameters |
| `a1`, `a2` | anisotropy scalings (unnormalized-Legendre convention) |
| `B_cm1` | rotational constant (cm⁻¹) |
| `M_amu` | molecular mass (amu) |
| `mu_He`, `mu_H2` | reduced masses for the two partners |
| `eps_0`, `eps_90`, `eps_180` | well depths at the three angles |
| `Rm_0`, `Rm_90`, `Rm_180` | minimum positions |
| `Rwall_0`, `Rwall_90`, `Rwall_180` | V=0 repulsive-wall radii |
| `d_eps_odd`, `d_eps_even` | anisotropy contrasts |
| `accepted` | (in a separate full log) bool + `reject_reason` |

- **Accepted library:** `potentials_library.csv`.
- **Full log (incl. rejects + reasons):** `potentials_log.csv` — for coverage
  analysis and range tuning.

---

## 9. Module / function layout (Python)

```
potential_generator.py
├── CONFIG                      # ranges, sample count, bounds, units, flags
├── v0(R, A, a, C6)             # Born–Mayer + dispersion
├── dv0_dR(R, A, a, C6)         # analytic derivative
├── angular_factor(costheta, a1, a2)   # unnormalized Legendre 1 + a1 P1 + a2 P2
├── find_minimum(A, a, C6, scale)      # outer root of dV/dR → (eps, Rm)
├── find_wall_radius(A, a, C6, scale)  # inner V=0 crossing
├── check_single_well(...)             # screening rule 1
├── check_no_collapse(...)             # screening rule 2
├── check_angular_physicality(a1, a2)  # screening rule 3 (full θ scan)
├── characterize(params)               # all descriptors at 0/90/180 + contrasts
├── sample_sobol(n)                    # scaled Sobol design over 7-D
├── generate(n)                        # loop: sample → screen → characterize
└── main()                             # write library CSV + full log CSV
```

### Dependencies
`numpy`, `scipy` (`stats.qmc.Sobol`, `optimize.brentq`), `pandas`.

---

## 10. Workflow position

```
[ this generator ]  →  potentials_library.csv
        │
        ▼
[ MOLSCAT ITYPE=1 input generator ]   (next component)
   - turns each library row into a MOLSCAT input deck
   - sets channel basis / energy grid scaled to each system's B and μ
        │
        ▼
[ scattering runs → σ(E) → k(T) fit: log k = Σ aₙ x^n, x = T^(−1/6) ]
        │
        ▼
[ regression / GP surrogate over potential space; active learning ]
```

---

## 11. Build order (suggested)

1. `v0`, `dv0_dR`, `angular_factor` — core math.
2. `find_minimum`, `find_wall_radius` — root finding + tests on a known case.
3. Screening functions + `characterize`.
4. `sample_sobol`, `generate`, CSV writers.
5. Run a small batch (e.g. 256), inspect coverage (§6), tune ranges.
6. Scale to the seed size (~200 accepted) for the active-learning start.
