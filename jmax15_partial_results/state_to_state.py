#!/usr/bin/env python3
"""State-to-state cross sections and rate coefficients (smoke test).

For each completed potential (He and H2):
  * full state-to-state inelastic cross sections sigma(i->f, E)  [Ang^2],
  * state-to-state rates k(i->f, T) via energy_kinteg_helper.k_of_T,
    integrating over the KINETIC energy in the initial channel,
        E_kin = E_total - E_i,   E_i = B*i*(i+1),
    (the shift only vanishes for i=0; it matters for every excited initial state,
    and makes de-excitation rates finite down to T->0),
  * fit  log10 k(i->f, T) = sum_n a_n x^n,  x = T^(-1/6).

Writes long tables (all transitions) + per-case wide tables from j=0.
"""

import os
import sys
import numpy as np
import pandas as pd

if not hasattr(np, "trapz"):          # NumPy 2.x compatibility for the helper
    np.trapz = np.trapezoid

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "smoke_test"))
import energy_kinteg_helper as ek          # noqa: E402
from parse_xs import parse_file            # noqa: E402

PROD = os.path.join(REPO, "production")
N_POT = 3
PARTNERS = ["He", "H2"]
T_GRID = np.geomspace(2.0, 50.0, 25)
FIT_DEG = 5
KMIN = 1.0e-12          # drop de-excitation transitions whose peak k(T) < this [cm^3/s]


def is_done(p):
    try:
        with open(p) as fh:
            return "This run used" in fh.read()
    except OSError:
        return False


def select_potentials(man):
    chosen = []
    for pid in man["id"].drop_duplicates():
        outs = {}
        ok = True
        for pn in PARTNERS:
            row = man[(man["id"] == pid) & (man["partner"] == pn)]
            if row.empty:
                ok = False
                break
            p = os.path.join(PROD, row.iloc[0]["dir"],
                             row.iloc[0]["input"].replace(".input", ".out"))
            if not is_done(p):
                ok = False
                break
            outs[pn] = (p, float(row.iloc[0]["ured"]), float(row.iloc[0]["B"]))
        if ok:
            chosen.append((pid, outs))
        if len(chosen) >= N_POT:
            break
    return chosen


def fit_logk(T, k):
    good = k > 0
    if good.sum() < FIT_DEG + 2:
        return [np.nan] * (FIT_DEG + 1), np.nan, int(good.sum())
    x = T[good] ** (-1.0 / 6.0)
    y = np.log10(k[good])
    coef = np.polyfit(x, y, FIT_DEG)
    a = coef[::-1]
    res = float(np.max(np.abs(y - np.polyval(coef, x))))
    return a, res, int(good.sum())


def main():
    man = pd.read_csv(os.path.join(PROD, "manifest.csv"))
    chosen = select_potentials(man)
    if not chosen:
        print("No potentials with both partners complete yet.")
        return
    print(f"State-to-state on {len(chosen)} potentials: "
          f"{[pid for pid, _ in chosen]}\n")

    xs_long, kt_long, coeff_rows = [], [], []
    n_kept = n_drop = 0
    for pid, outs in chosen:
        for pn, (out, mu, B) in outs.items():
            df = parse_file(out)
            # de-excitation only: final j below initial j (exothermic, the useful ones)
            deexc = df[df["j_F"] < df["j_I"]].copy()

            for (i, f), sub in deexc.groupby(["j_I", "j_F"]):
                sub = sub.sort_values("E")
                Ekin = sub["E"].values - B * i * (i + 1)   # kinetic E in channel i
                m = Ekin > 0
                if m.sum() < 3:
                    continue
                Ek, sg = Ekin[m], sub["SIG"].values[m]
                k = np.array([ek.k_of_T(Ek, sg, T, mu) for T in T_GRID])
                if not np.isfinite(k).any() or np.nanmax(k) < KMIN:
                    n_drop += 1                            # negligible transition
                    continue
                n_kept += 1

                # cross sections for this (kept) transition
                xs = (sub.loc[m, ["E", "j_I", "j_F", "SIG"]]
                      .rename(columns={"E": "E_cm1", "j_I": "j_i",
                                       "j_F": "j_f", "SIG": "sigma_ang2"}))
                xs.insert(0, "partner", pn)
                xs.insert(0, "id", pid)
                xs_long.append(xs)

                kt_long.append(pd.DataFrame({"id": pid, "partner": pn,
                                             "j_i": i, "j_f": f,
                                             "T_K": T_GRID, "k_cm3_s": k}))
                a, res, npts = fit_logk(T_GRID, k)
                row = {"id": pid, "partner": pn, "j_i": i, "j_f": f,
                       "mu_amu": round(mu, 4), "k_peak": float(np.nanmax(k)),
                       "n_fit_pts": npts,
                       "max_resid_log10": (round(res, 4) if res == res else np.nan)}
                row.update({f"a{n}": a[n] for n in range(len(a))})
                coeff_rows.append(row)

    pd.concat(xs_long).to_csv(os.path.join(HERE, "s2s_xs_all.csv"), index=False)
    pd.concat(kt_long).to_csv(os.path.join(HERE, "s2s_kT_all.csv"), index=False)
    coeffs = pd.DataFrame(coeff_rows)
    coeffs.to_csv(os.path.join(HERE, "s2s_kT_coeffs.csv"), index=False)

    # ---- console report -------------------------------------------------- #
    print(f"de-excitation transitions kept (peak k >= {KMIN:.0e}): {n_kept}   "
          f"dropped (too small): {n_drop}   "
          f"cross-section rows: {len(pd.concat(xs_long))}")

    ex = chosen[0][0]
    sub = coeffs[(coeffs["id"] == ex) & (coeffs["partner"] == "He")]
    print(f"\nExample pot{ex:04d} He - de-excitation k(i->f) [cm^3/s]:")
    kt = pd.concat(kt_long)
    kt = kt[(kt["id"] == ex) & (kt["partner"] == "He")]
    Tsamp = [5, 20, 50]
    pairs = [(1, 0), (2, 0), (2, 1), (3, 1), (3, 2)]
    print("  i->f " + "".join(f"{t:>14}K" for t in Tsamp))
    for i, f in pairs:
        s = kt[(kt["j_i"] == i) & (kt["j_f"] == f)]
        if s.empty:
            continue
        vals = "".join(f"{np.interp(t, s['T_K'], s['k_cm3_s']):>15.4e}"
                       for t in Tsamp)
        print(f"  {i}->{f}{vals}")

    print(f"\nFit coeffs (log10 k = sum a_n x^n, x=T^-1/6) sample, pot{ex:04d} He:")
    cols = ["j_i", "j_f", "k_peak", "n_fit_pts", "max_resid_log10"] + \
           [f"a{n}" for n in range(FIT_DEG + 1)]
    show = sub[(sub["j_i"] <= 3)].sort_values(["j_i", "j_f"])
    with pd.option_context("display.width", 220, "display.max_columns", 25,
                           "display.float_format", lambda v: f"{v:.3g}"):
        print(show[cols].to_string(index=False))

    print(f"\nWrote s2s_xs_all.csv, s2s_kT_all.csv, s2s_kT_coeffs.csv to {HERE}")


if __name__ == "__main__":
    main()
