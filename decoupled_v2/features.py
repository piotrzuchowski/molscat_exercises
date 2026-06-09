#!/usr/bin/env python3
"""Physics-informed coupling descriptors for the potential library.

For each row we summarize the anisotropic radial functions V1, V2 (the things
that actually drive |dj|=1, 2 transitions) by integral transforms over the
dynamically accessible region [Rwall, R_out] -- NOT pointwise grid values.
These are cheap, differentiable, low-dimensional, and meant to *correlate* with
inelastic cross sections (they need not predict them).

Per lambda in {1,2}:
    I0_l = \int V_l dR                       zeroth moment (net coupling)
    I1_l = \int V_l R dR                     first moment (where it acts)
    Iw_l = \int V_l / kappa(R) dR            semiclassical dwell-time weight
    Fq_l = \int V_l j_l(qR) R^2 dR           Born form factor at thermal q

kappa and q use a thermal reference energy E_REF ~ k_B*50 K.  Units: V in cm^-1,
R in Ang.  Run:  python features.py  ->  features.csv
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import spherical_jn

HARTREE2CM = 219474.6313702
BOHR2ANG = 0.529177210903
HBAR2_2MU = 16.85762920

E_REF = 35.0       # cm^-1  (~ k_B * 50 K) reference collision energy
R_OUT = 15.0       # Ang outer integration limit
N_R = 3000
MU = 4.002602 * 20.0 / (4.002602 + 20.0)   # He + mass-20 molecule [amu] = 3.335


def v_lambda(R, A_au, a_au, C6_au, a1, a2, eta):
    """Return (V0, V1, V2) in cm^-1 on grid R [Ang]."""
    A = A_au * HARTREE2CM
    a = a_au / BOHR2ANG
    C6 = C6_au * HARTREE2CM * BOHR2ANG ** 6
    exp = A * np.exp(-a * R)
    disp = C6 * R ** -6
    return exp - disp, eta * exp - a1 * disp, a2 * (exp - disp)


def features_row(row):
    Rin = float(row["Rwall_90"])                 # start where the wave turns on
    R = np.linspace(Rin, R_OUT, N_R)
    V0, V1, V2 = v_lambda(R, row["A"], row["a"], row["C6"],
                          row["a1"], row["a2"], row["eta"])

    k2 = E_REF * MU / HBAR2_2MU                   # thermal wavenumber^2 [Ang^-2]
    kappa = np.sqrt(np.clip(k2 - V0 * MU / HBAR2_2MU, 0.1 ** 2, None))
    q = 2.0 * np.sqrt(k2)                         # backscattering momentum transfer

    out = {"id": int(row["id"])}
    for lam, Vl in ((1, V1), (2, V2)):
        out[f"I0_{lam}"] = trapezoid(Vl, R)
        out[f"I1_{lam}"] = trapezoid(Vl * R, R)
        out[f"Iw_{lam}"] = trapezoid(Vl / kappa, R)
        out[f"Fq_{lam}"] = trapezoid(Vl * spherical_jn(lam, q * R) * R ** 2, R)
    return out


def main(lib_csv="potentials_library.csv", out_csv="features.csv"):
    df = pd.read_csv(lib_csv)
    feats = pd.DataFrame([features_row(r) for _, r in df.iterrows()])
    feats.to_csv(out_csv, index=False)
    print(f"Wrote {len(feats)} rows x {feats.shape[1]-1} features -> {out_csv}")
    print(feats.describe().loc[["mean", "std", "min", "max"]].T.to_string())
    return feats


if __name__ == "__main__":
    main()
