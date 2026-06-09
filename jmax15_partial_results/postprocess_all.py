#!/usr/bin/env python3
"""Post-process every completed potential: de-excitation state-to-state cross
sections + k(T) fit coefficients. Streams to disk (append per potential) so a
crash/kill never loses finished work. Reuses the validated state_to_state core.

Outputs (overwritten fresh each run):
  s2s_xs_all.csv      id, partner, E_cm1, j_i, j_f, sigma_ang2   (de-excitation)
  s2s_kT_coeffs.csv   id, partner, j_i, j_f, mu_amu, k_peak, n_fit_pts,
                      max_resid_log10, a0..a5   (log10 k = sum a_n (T^-1/6)^n)
"""

import os
import sys
import numpy as np
import pandas as pd

if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "smoke_test"))
import energy_kinteg_helper as ek          # noqa: E402
from parse_xs import parse_file            # noqa: E402
from state_to_state import (fit_logk, is_done, T_GRID, FIT_DEG, KMIN)  # noqa: E402

PROD = os.path.join(REPO, "production")
XS_CSV = os.path.join(HERE, "s2s_xs_all.csv")
COEF_CSV = os.path.join(HERE, "s2s_kT_coeffs.csv")
PARTNERS = ["He", "H2"]


def complete_potentials(man):
    out = []
    for pid in man["id"].drop_duplicates():
        rec = {}
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
            rec[pn] = (p, float(row.iloc[0]["ured"]), float(row.iloc[0]["B"]))
        if ok:
            out.append((int(pid), rec))
    return out


def process_one(pid, pn, out, mu, B):
    """Return (xs_df, coeff_rows) for one deck's de-excitation transitions."""
    df = parse_file(out)
    deexc = df[df["j_F"] < df["j_I"]]
    xs_parts, coeffs = [], []
    for (i, f), sub in deexc.groupby(["j_I", "j_F"]):
        sub = sub.sort_values("E")
        Ekin = sub["E"].values - B * i * (i + 1)
        m = Ekin > 0
        if m.sum() < 3:
            continue
        Ek, sg = Ekin[m], sub["SIG"].values[m]
        k = np.array([ek.k_of_T(Ek, sg, T, mu) for T in T_GRID])
        if not np.isfinite(k).any() or np.nanmax(k) < KMIN:
            continue
        xs = (sub.loc[m, ["E", "j_I", "j_F", "SIG"]]
              .rename(columns={"E": "E_cm1", "j_I": "j_i",
                               "j_F": "j_f", "SIG": "sigma_ang2"}))
        xs.insert(0, "partner", pn)
        xs.insert(0, "id", pid)
        xs_parts.append(xs)
        a, res, npts = fit_logk(T_GRID, k)
        row = {"id": pid, "partner": pn, "j_i": int(i), "j_f": int(f),
               "mu_amu": round(mu, 5), "k_peak": float(np.nanmax(k)),
               "n_fit_pts": npts,
               "max_resid_log10": (round(res, 4) if res == res else np.nan)}
        row.update({f"a{n}": a[n] for n in range(FIT_DEG + 1)})
        coeffs.append(row)
    xs_df = pd.concat(xs_parts) if xs_parts else pd.DataFrame()
    return xs_df, coeffs


def main():
    man = pd.read_csv(os.path.join(PROD, "manifest.csv"))
    pots = complete_potentials(man)
    print(f"Processing {len(pots)} complete potentials "
          f"(both partners), KMIN={KMIN:.0e}", flush=True)
    for f in (XS_CSV, COEF_CSV):
        if os.path.exists(f):
            os.remove(f)

    n_pot = n_trans = 0
    xs_hdr = coef_hdr = False
    for pid, rec in pots:
        xs_acc, coef_acc = [], []
        for pn in PARTNERS:
            out, mu, B = rec[pn]
            xs_df, coeffs = process_one(pid, pn, out, mu, B)
            if not xs_df.empty:
                xs_acc.append(xs_df)
            coef_acc.extend(coeffs)
            n_trans += len(coeffs)
        if xs_acc:
            pd.concat(xs_acc).to_csv(XS_CSV, mode="a", header=not xs_hdr,
                                     index=False)
            xs_hdr = True
        if coef_acc:
            pd.DataFrame(coef_acc).to_csv(COEF_CSV, mode="a", header=not coef_hdr,
                                          index=False)
            coef_hdr = True
        n_pot += 1
        if n_pot % 25 == 0:
            print(f"  {n_pot}/{len(pots)} potentials, "
                  f"{n_trans} transitions so far", flush=True)

    print(f"DONE: {n_pot} potentials, {n_trans} de-excitation transitions kept.",
          flush=True)
    print(f"  -> {XS_CSV}")
    print(f"  -> {COEF_CSV}")


if __name__ == "__main__":
    main()
