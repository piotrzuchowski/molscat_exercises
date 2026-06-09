#!/usr/bin/env python3
"""Smoke test of the analysis pipeline on completed production decks.

For a few potentials whose He AND H2 runs are finished:
  1. parse state-to-state integral cross sections,
  2. build the total inelastic (rotational-excitation) cross section from j=0,
     sigma_inel(E) = sum_{f!=0} sigma(0->f),
  3. thermally average -> rate k(T) via energy_kinteg_helper.k_of_T,
  4. fit  log10 k(T) = sum_n a_n * x^n   with x = T^(-1/6),
and writes tables + coefficients into this directory.
"""

import os
import sys
import glob
import numpy as np
import pandas as pd

# NumPy 2.x removed np.trapz (renamed np.trapezoid); energy_kinteg_helper uses
# np.trapz. Restore the alias on the shared module so the helper works unchanged.
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "smoke_test"))
import energy_kinteg_helper as ek          # noqa: E402
from parse_xs import parse_file            # noqa: E402

PROD = os.path.join(REPO, "production")
N_POT = 3                                  # potentials to smoke-test
PARTNERS = ["He", "H2"]
T_GRID = np.geomspace(2.0, 50.0, 25)       # K (helper grid supports up to ~50 K)
FIT_DEG = 5                                # log10 k = sum_{0..DEG} a_n x^n


def is_done(path):
    try:
        with open(path) as fh:
            return "This run used" in fh.read()
    except OSError:
        return False


def sigma_inel_from_ground(df):
    """Full-grid (E, sigma_inel) with sigma=0 below the first open threshold."""
    E = np.sort(df["E"].unique())
    g = df[(df["I"] == 1) & (df["F"] != 1)]
    s = g.groupby("E")["SIG"].sum().reindex(E).fillna(0.0).values
    return E, s


def fit_logk(T, k):
    """Fit log10 k = sum a_n x^n, x=T^-1/6. Returns a[0..DEG] (a0=constant)."""
    good = k > 0
    x = T[good] ** (-1.0 / 6.0)
    y = np.log10(k[good])
    coef = np.polyfit(x, y, FIT_DEG)        # highest power first
    a = coef[::-1]                          # a0..aDEG
    resid = y - np.polyval(coef, x)
    return a, float(np.max(np.abs(resid))), good.sum()


def select_potentials(man):
    """First N_POT library ids whose He and H2 outputs are both complete."""
    chosen = []
    for pid in man["id"].drop_duplicates():
        outs = {pn: None for pn in PARTNERS}
        for pn in PARTNERS:
            row = man[(man["id"] == pid) & (man["partner"] == pn)]
            if row.empty:
                break
            p = os.path.join(PROD, row.iloc[0]["dir"],
                             row.iloc[0]["input"].replace(".input", ".out"))
            if is_done(p):
                outs[pn] = p
        if all(outs.values()):
            chosen.append((pid, outs))
        if len(chosen) >= N_POT:
            break
    return chosen


def main():
    man = pd.read_csv(os.path.join(PROD, "manifest.csv"))
    chosen = select_potentials(man)
    if not chosen:
        print("No potentials with both partners complete yet.")
        return
    print(f"Smoke test on {len(chosen)} potentials: "
          f"{[pid for pid, _ in chosen]}\n")

    inel_tables, coeff_rows, kt_long = [], [], []
    for pid, outs in chosen:
        for pn in PARTNERS:
            mu = float(man[(man["id"] == pid) &
                           (man["partner"] == pn)].iloc[0]["ured"])
            df = parse_file(outs[pn])
            E, s = sigma_inel_from_ground(df)

            # inelastic cross section table
            it = pd.DataFrame({"E_cm1": E, "sigma_inel_ang2": s})
            it.insert(0, "partner", pn)
            it.insert(0, "id", pid)
            it.to_csv(os.path.join(HERE, f"inelastic_pot{pid:04d}_{pn}.csv"),
                      index=False)
            inel_tables.append(it)

            # rate k(T)
            k = np.array([ek.k_of_T(E, s, T, mu) for T in T_GRID])
            kt = pd.DataFrame({"T_K": T_GRID, "k_cm3_s": k})
            kt.insert(0, "partner", pn)
            kt.insert(0, "id", pid)
            kt.to_csv(os.path.join(HERE, f"kT_pot{pid:04d}_{pn}.csv"),
                      index=False)
            kt_long.append(kt)

            a, maxres, npts = fit_logk(T_GRID, k)
            row = {"id": pid, "partner": pn, "mu_amu": round(mu, 4),
                   "n_fit_pts": npts, "max_resid_log10": round(maxres, 4)}
            row.update({f"a{n}": a[n] for n in range(len(a))})
            coeff_rows.append(row)

    coeffs = pd.DataFrame(coeff_rows)
    coeffs.to_csv(os.path.join(HERE, "kT_fit_coeffs.csv"), index=False)
    pd.concat(inel_tables).to_csv(os.path.join(HERE, "inelastic_all.csv"),
                                  index=False)
    pd.concat(kt_long).to_csv(os.path.join(HERE, "kT_all.csv"), index=False)

    # ---- console report -------------------------------------------------- #
    print("Inelastic cross section sigma_inel(0->all) [Ang^2] at sample E:")
    Esamp = [1, 10, 50, 100, 200]
    hdr = "  id    partner " + "".join(f"{e:>10}cm" for e in Esamp)
    print(hdr)
    for it in inel_tables:
        E = it["E_cm1"].values
        s = it["sigma_inel_ang2"].values
        vals = "".join(f"{np.interp(e, E, s):>12.4g}" for e in Esamp)
        print(f"  {it['id'].iloc[0]:<5} {it['partner'].iloc[0]:<6}{vals}")

    print("\nk(T) [cm^3/s] at sample T:")
    Tsamp = [5, 10, 20, 50]
    print("  id    partner " + "".join(f"{t:>13}K" for t in Tsamp))
    for kt in kt_long:
        vals = "".join(f"{np.interp(t, kt['T_K'], kt['k_cm3_s']):>14.4e}"
                       for t in Tsamp)
        print(f"  {kt['id'].iloc[0]:<5} {kt['partner'].iloc[0]:<6}{vals}")

    print(f"\nFit  log10 k(T) = sum_n a_n * x^n,  x = T^(-1/6),  deg={FIT_DEG}:")
    with pd.option_context("display.width", 200, "display.max_columns", 20,
                           "display.float_format", lambda v: f"{v:.4g}"):
        print(coeffs.to_string(index=False))
    print(f"\nWrote tables + kT_fit_coeffs.csv to {HERE}")


if __name__ == "__main__":
    main()
