# decoupled_v3 — V1 sampled by MAGNITUDE (odd anisotropy strength)

## Why this version
Under reversal of the molecular axis, `theta -> pi - theta`, the Legendre
polynomials transform as `P_lambda -> (-1)^lambda P_lambda`. Hence

    V1 -> -V1 ,   V2 -> V2 .

The **sign of V1 is not physical** — it only reflects the arbitrary choice of
which end of the rotor is `theta = 0`. For collision rates of unoriented linear
rotors the meaningful quantity is the **strength** of the odd anisotropy,
`|V1|/|V0|`, not its sign. Two potentials differing only by `V1 -> -V1` give
identical rates and must give identical ML descriptors.

decoupled_v3 therefore (i) samples the *magnitude* `rho1 = |V1|/|V0| in [0,1]`,
(ii) gives `V1` its own independent radial form, and (iii) uses absolute-value
odd-anisotropy descriptors. `V2` (even, sign IS physical) is unchanged from v2.

## Potential model
Unnormalised Legendre expansion (`P1 = cos t`, `P2 = 3cos^2 t - 1`):

    V(R,theta) = V0(R) + s1 * |V1|(R) * P1 + V2(R) * P2

    V0(R) =      A e^{-aR}      -      C6 R^{-6}         (isotropic)
    V2(R) = zeta A e^{-aR}      - a2   C6 R^{-6}         (even; decoupled wall/tail)
    |V1|(R) = A1 e^{-b1 R} - C7 R^{-7}                   (odd MAGNITUDE, R^-7 tail)

- `A, a, C6` in atomic units (Hartree, 1/Bohr, Hartree·Bohr^6).
- The odd term has the **dipole-type R^-7** long-range dispersion (not R^-6), its
  **own steepness `b1`** (independent of `a`), and is multiplied by an internal
  sign `s1 = ±1`.
- `s1` is a *geometry-only* gauge: it is carried into the MOLSCAT deck so the
  drawn potential is concrete, but it is **never** an ML feature.

## How V1 is parametrized
`|V1|(R)` is built from *physical* parameters so it has a well of depth `eps1`
at position `R_{m,1}`, with reduced steepness `x1 = b1 R_{m,1}`:

    eps1 = rho1 * eps0                 (eps0 = isotropic V0 well depth)
    R_{m,1} = m1 * R_m(V0)             (V1 well placed near the V0 well)
    b1 = x1 / R_{m,1}
    A1 = 7 eps1 e^{x1} / (x1 - 7)
    C7 = eps1 x1 R_{m,1}^7 / (x1 - 7)

These solve `|V1|(R_{m,1}) = -eps1`, `|V1|'(R_{m,1}) = 0`. Positivity of `A1`
requires `x1 > 7`; we sample `x1 ∈ [9,15]` (default; stress-test up to 18).

### Sampled parameters (10-D Sobol + 1 sign)
| symbol | meaning | range / law |
|---|---|---|
| `A` | repulsive prefactor [Hartree] | log-uniform [150, 5000] |
| `a` | repulsive steepness [1/Bohr] | log-uniform [2.2, 2.9] |
| `C6` | isotropic dispersion [au] | log-uniform [3, 40] |
| `a2` | even tail anisotropy | uniform [−0.4, 0.4] |
| `zeta` | even wall anisotropy | Gaussian, σ=0.2 |
| `B` | rotational constant [cm⁻¹] | log-uniform [1e-3, 10] |
| `M` | mass [amu] | log-uniform [8, 100] |
| **`rho1`** | **odd-anisotropy strength `\|V1\|/\|V0\|`** | **mixture (below)** |
| `m1` | `R_{m,1}/R_m(V0)` | uniform [0.90, 1.15] |
| `x1` | reduced steepness `b1 R_{m,1}` | uniform [9, 15] |
| `s1` | internal odd sign | ±1 (rng, not a feature) |

### rho1 mixture — oversample weak odd anisotropy
Weak `V1` produces small, hard-to-learn Δj=1 rates (esp. 1→0), so we oversample
the low end:

- 40% : `rho1 ∈ [1e-3, 0.05]`, log-uniform
- 40% : `rho1 ∈ [0.05, 0.5]`,  log-uniform
- 20% : `rho1 ∈ [0.5, 1.0]`,   uniform

## ML descriptors (absolute-valued)
Replace the signed odd contrast by magnitudes (`eps_floor = 1e-6 cm⁻¹`):

    abs_d_eps_odd = |eps_0 - eps_180|
    rel_odd       = |eps_0 - eps_180| / (eps_90 + eps_floor)

The signed `d_eps_odd_signed` is kept **only for diagnostics** ("which end is
more attractive") and must not be used as a learning feature. Even-anisotropy
descriptor `d_eps_even = (eps_0+eps_180)/2 - eps_90` is unchanged.

**Consistency check:** two potentials differing only by `V1 -> -V1` produce
identical `abs_d_eps_odd`, `rel_odd` (and identical rates up to numerical noise).

## MOLSCAT deck mapping
The λ=1 term becomes `NPOWER = 0, -7` with its own exponent:

    LAMBDA  = 0,  1,  2
    NPOWER  = 0, -6,   0, -7,   0, -6
    E (exp) = -a,      -b1,     -a            (per NPOWER==0 term)
    A       = A0, -C6_0,  s1*A1, -s1*C7,  zeta*A0, -a2*C6_0   (cm⁻¹, Å)

Units: `A1[Hartree]·HARTREE2CM`, `b1[1/Bohr]/BOHR2ANG`,
`C7[Hartree·Bohr^7]·HARTREE2CM·BOHR2ANG^7`.

## Validation (smoke test)
σ(1→0) tracks `rho1` monotonically over 3+ decades (≈0.016→36 Å² as
rho1: 0.01→0.88), while σ(2→0) stays decoupled (driven by `zeta`/`a2`).
Confirms odd-anisotropy strength controls Δj=1 and is independent of Δj=2.

## Campaign
- Library: `potential_generator.py` (seed=2026, sobol_m=15) → 15,200 accepted
  (≈46% acceptance). Columns include `rho1, m1, x1, s1, A1_au, b1_au, C7_au,
  Rm1_ang, eps1_cm1, abs_d_eps_odd, rel_odd`.
- Run: first **10,000** potentials, H2 partner, **jmax=15**, dense k(T) grid
  (0.05–208.5 cm⁻¹, jmax15 machinery), DR=0.04, DTOL=10, BCYOMN hard wall.
- Build/submit: `build_campaign.py [n]` → `production/`; `run_all.sh`
  (offset-chunked ≤180 concurrent); `postprocess.py` →
  `molscat_xsec.csv` + `s2s_kT_coeffs.csv`.
