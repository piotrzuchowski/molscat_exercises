#!/usr/bin/env python3
"""Self-contained evaluator for the model potentials in potentials_library.csv.

Potential model (atom + linear rigid rotor), expanded in UNNORMALISED Legendre
polynomials  P1 = cos(theta),  P2 = 3 cos^2(theta) - 1 :

    V(R, theta) = V0(R) + V1(R) P1(cos t) + V2(R) P2(cos t)

    V0(R) =      A e^{-aR} -      C6 R^{-6}        (isotropic)
    V1(R) =  eta A e^{-aR} - a1   C6 R^{-6}        (lambda = 1 anisotropy)
    V2(R) = a2 ( A e^{-aR} -      C6 R^{-6}) = a2 * V0(R)   (lambda = 2)

UNITS.  The columns A, a, C6 in potentials_library.csv are ATOMIC UNITS:
    A  in Hartree,  a in 1/Bohr,  C6 in Hartree * Bohr^6.
So V0/V1/V2 below take R in BOHR and return energy in HARTREE.  Use the
constants to convert:  E[cm^-1] = E[Hartree] * HARTREE2CM ,
                       R[Ang]  = R[Bohr]    * BOHR2ANG .

The DWBA-style feature idea (what this file is for): a lambda=1 / lambda=2
inelastic transition amplitude in distorted-wave Born approximation scales with
the radial overlap of V_lambda(R) with the incoming/outgoing distorted waves,
which live mostly inside the V0 well.  A cheap, basis-set-free proxy is the
integral of V1(R) and V2(R) over the V0 well.  `well_bounds` + `dwba_features`
below give one concrete realisation; tweak the weight/limits as you like.
"""
import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import quad

HARTREE2CM = 219474.6313702      # Hartree -> cm^-1
BOHR2ANG = 0.529177210903        # Bohr -> Angstrom


# --------------------------------------------------------------------------- #
# Legendre components (R in Bohr, energy in Hartree)                           #
# --------------------------------------------------------------------------- #
def v0(R, A, a, C6):
    return A * np.exp(-a * R) - C6 * R ** -6


def v1(R, A, a, C6, a1, eta):
    return eta * A * np.exp(-a * R) - a1 * C6 * R ** -6


def v2(R, A, a, C6, a2):
    return a2 * (A * np.exp(-a * R) - C6 * R ** -6)


def v_theta(R, A, a, C6, a1, a2, eta, costheta):
    """Full V(R, theta); costheta in [-1, 1]."""
    p1 = costheta
    p2 = 3.0 * costheta ** 2 - 1.0
    return (A * np.exp(-a * R) * (1.0 + eta * p1 + a2 * p2)
            - C6 * R ** -6 * (1.0 + a1 * p1 + a2 * p2))


# --------------------------------------------------------------------------- #
# V0 well geometry (all radii returned in Bohr, depths in cm^-1)              #
# --------------------------------------------------------------------------- #
def well_bounds(A, a, C6):
    """Characterise the isotropic V0 well.

    Returns dict with (Bohr): R_in  (inner zero crossing V0=0, repulsive side),
    R_min (minimum), R_out is +inf (V0->0^- from below as R->inf).
    Also eps0_cm1 (well depth, positive), R_min_ang.
    The 'well' for integration is taken as R in [R_in, R_cut] where the
    integrand has decayed; callers pass their own R_cut (default 30 Bohr).
    Returns None if there is no bound well in the searched range.
    """
    # minimum of V0 between the repulsive wall and long range
    res = minimize_scalar(lambda R: v0(R, A, a, C6), bounds=(2.0, 40.0),
                          method="bounded")
    R_min = float(res.x)
    if v0(R_min, A, a, C6) >= 0:
        return None                       # purely repulsive, no well
    # inner zero crossing: V0 changes sign between wall and R_min
    try:
        R_in = brentq(v0, 1.0, R_min, args=(A, a, C6))
    except ValueError:
        R_in = 1.0
    eps0 = -v0(R_min, A, a, C6) * HARTREE2CM
    return {"R_in": R_in, "R_min": R_min, "eps0_cm1": eps0,
            "R_min_ang": R_min * BOHR2ANG}


# --------------------------------------------------------------------------- #
# Example DWBA-spirit features (edit freely)                                   #
# --------------------------------------------------------------------------- #
def dwba_features(row, R_cut=30.0, weight="depth"):
    """Integrals of V1, V2 over the V0 well -> a few scalar features.

    row : a pandas Series / dict with A, a, C6, a1, a2, eta.
    R_cut : outer integration limit (Bohr).
    weight :
        'none'  -> plain  integral_{R_in}^{R_cut} V_lambda(R) dR
        'depth' -> overlap with the V0 well shape:
                   integral V_lambda(R) * [-V0(R)]_+ dR   (DWBA-like overlap)
        'R2'    -> volume element:  integral V_lambda(R) R^2 dR

    All integrals are done in atomic units (R Bohr, V Hartree) then the energy
    factor is converted to cm^-1; lengths stay in Bohr (so weight='none' units
    are cm^-1 * Bohr, 'R2' are cm^-1 * Bohr^3, 'depth' are cm^-1^2 * Bohr).
    Returns dict of scalar features (NaN if no well).
    """
    A, a, C6 = row["A"], row["a"], row["C6"]
    a1, a2, eta = row["a1"], row["a2"], row["eta"]
    wb = well_bounds(A, a, C6)
    if wb is None:
        return {k: np.nan for k in
                ("I_V1", "I_V2", "I_V0", "I_V1_norm", "I_V2_norm")}
    R_in = wb["R_in"]

    def wfac(R):
        if weight == "none":
            return 1.0
        if weight == "R2":
            return R ** 2
        # 'depth': weight by the (positive) well depth profile
        v = v0(R, A, a, C6)
        return max(-v, 0.0)

    def integ(fun):
        val, _ = quad(lambda R: fun(R) * wfac(R) * HARTREE2CM,
                      R_in, R_cut, limit=200)
        return val

    I_V0 = integ(lambda R: v0(R, A, a, C6))
    I_V1 = integ(lambda R: v1(R, A, a, C6, a1, eta))
    I_V2 = integ(lambda R: v2(R, A, a, C6, a2))
    # normalised by the V0 integral -> dimensionless anisotropy of the overlap
    norm = I_V0 if I_V0 != 0 else np.nan
    return {"I_V0": I_V0, "I_V1": I_V1, "I_V2": I_V2,
            "I_V1_norm": I_V1 / norm, "I_V2_norm": I_V2 / norm}


def build_feature_table(lib_csv="potentials_library.csv", **kw):
    """Compute DWBA features for every potential -> DataFrame keyed by id."""
    lib = pd.read_csv(lib_csv)
    rows = []
    for _, r in lib.iterrows():
        f = dwba_features(r, **kw)
        f["id"] = int(r["id"])
        rows.append(f)
    return pd.DataFrame(rows).set_index("id")


if __name__ == "__main__":
    # demo: print well + DWBA features for the first few potentials
    lib = pd.read_csv("potentials_library.csv")
    print("id   eps0[cm-1]  Rmin[Ang]   I_V1(depth)   I_V2(depth)  "
          "I_V1/V0   I_V2/V0")
    for _, r in lib.head(8).iterrows():
        wb = well_bounds(r["A"], r["a"], r["C6"])
        f = dwba_features(r)
        if wb is None:
            print(f"{int(r['id']):>3}  (no well)")
            continue
        print(f"{int(r['id']):>3}  {wb['eps0_cm1']:>9.2f}  "
              f"{wb['R_min_ang']:>9.3f}  {f['I_V1']:>12.3e}  "
              f"{f['I_V2']:>12.3e}  {f['I_V1_norm']:>8.3f}  "
              f"{f['I_V2_norm']:>8.3f}")
