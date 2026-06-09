#!/usr/bin/env python3
"""Parse every finished decoupled_v2 deck -> cross sections + k(T) fits.

Outputs:
  molscat_xsec.csv    id, B_cm1, mu, energy, I, F, J_I, J_F, kind, sigma
  s2s_kT_coeffs.csv   id, partner, j_i, j_f, mu_amu, k_peak, n_fit_pts,
                      max_resid_log10, a0..a5   (log10 k = sum a_n (T^-1/6)^n)
  run_status.csv      id, status, n_sig

Reuses the jmax15 k(T) machinery: de-excitation transitions (j_f<j_i), kinetic
energy E_kin = E - B*i(i+1), Boltzmann average over T_GRID, T^-1/6 polynomial.
"""
import os
import sys
import glob
import numpy as np
import pandas as pd

if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

HERE = os.path.dirname(os.path.abspath(__file__))
DV2 = os.path.dirname(HERE)
REPO = os.path.dirname(DV2)
sys.path.insert(0, DV2)
sys.path.insert(0, REPO)
import molscat_run as mr
import energy_kinteg_helper as ek

PARTNER = "H2"                       # campaign uses mu_H2
T_GRID = np.geomspace(2.0, 50.0, 25)
FIT_DEG = 5
KMIN = 1.0e-12
XSEC = os.path.join(HERE, "molscat_xsec.csv")
COEF = os.path.join(HERE, "s2s_kT_coeffs.csv")
STATUS = os.path.join(HERE, "run_status.csv")


def fit_logk(T, k):
    good = k > 0
    if good.sum() < FIT_DEG + 2:
        return [np.nan] * (FIT_DEG + 1), np.nan, int(good.sum())
    x = T[good] ** (-1.0 / 6.0)
    y = np.log10(k[good])
    coef = np.polyfit(x, y, FIT_DEG)
    resid = float(np.max(np.abs(np.polyval(coef, x) - y)))
    return list(coef[::-1]), resid, int(good.sum())   # a0..a5


def kT_rows(sig_df, cid, B, mu):
    """De-excitation k(T) fits for one potential."""
    out = []
    deexc = sig_df[sig_df["J_F"] < sig_df["J_I"]]
    for (i, f), sub in deexc.groupby(["J_I", "J_F"]):
        sub = sub.sort_values("energy")
        Ekin = sub["energy"].values - B * i * (i + 1)
        m = Ekin > 0
        if m.sum() < 3:
            continue
        Ek, sg = Ekin[m], sub["sigma"].values[m]
        k = np.array([ek.k_of_T(Ek, sg, T, mu) for T in T_GRID])
        if not np.isfinite(k).any() or np.nanmax(k) < KMIN:
            continue
        a, resid, npts = fit_logk(T_GRID, k)
        row = {"id": cid, "partner": PARTNER, "j_i": int(i), "j_f": int(f),
               "mu_amu": round(mu, 5), "k_peak": float(np.nanmax(k)),
               "n_fit_pts": npts,
               "max_resid_log10": (round(resid, 4) if resid == resid else np.nan)}
        row.update({f"a{n}": a[n] for n in range(FIT_DEG + 1)})
        out.append(row)
    return out


def main():
    man = pd.read_csv(os.path.join(HERE, "manifest.csv")).set_index("id")
    xs_rows, coef_rows, status = [], [], []
    for out in sorted(glob.glob(os.path.join(HERE, "decks", "pot_*.out"))):
        cid = int(os.path.basename(out)[4:9])
        txt = open(out).read()
        if "This run used" not in txt:
            status.append({"id": cid, "status": "incomplete", "n_sig": 0})
            continue
        sig = mr.parse_sigma(txt)
        if not sig:
            status.append({"id": cid, "status": "no_xsec", "n_sig": 0})
            continue
        B = float(man.loc[cid, "B_cm1"]); mu = float(man.loc[cid, "mu"])
        df = pd.DataFrame(sig)
        df["id"] = cid; df["B_cm1"] = B; df["mu"] = mu
        df["kind"] = np.where(df["F"] == df["I"], "elastic", "inelastic")
        xs_rows.append(df[["id", "B_cm1", "mu", "energy", "I", "F",
                           "J_I", "J_F", "kind", "sigma"]])
        coef_rows.extend(kT_rows(df, cid, B, mu))
        status.append({"id": cid, "status": "ok", "n_sig": len(sig)})

    if xs_rows:
        xs = pd.concat(xs_rows, ignore_index=True)
        xs.to_csv(XSEC, index=False)
    if coef_rows:
        pd.DataFrame(coef_rows).to_csv(COEF, index=False)
    st = pd.DataFrame(status)
    st.to_csv(STATUS, index=False)

    nok = (st["status"] == "ok").sum()
    print(f"potentials: {len(st)}  ok={nok}  "
          f"incomplete={(st['status']=='incomplete').sum()}  "
          f"no_xsec={(st['status']=='no_xsec').sum()}")
    if xs_rows:
        print(f"cross sections: {len(xs)} "
              f"(inelastic={(xs['kind']=='inelastic').sum()}) -> {XSEC}")
    if coef_rows:
        print(f"k(T) fits: {len(coef_rows)} de-excitation transitions -> {COEF}")
    print(f"status -> {STATUS}")


if __name__ == "__main__":
    main()
