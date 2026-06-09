#!/usr/bin/env python3
"""Fingerprint vs physics-derived features, across the dj=1 ladder + H2 partner.

Part A: target He 1->0 -- compare the raw V(R,theta) fingerprint against derived
        geometry and DWBA-overlap feature sets (GBT and GP), same 5-fold CV.
Part B: sweep (He 1->0, 2->1, 3->2; H2 1->0) with GBT to see where DWBA helps.
"""
import warnings
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
from sklearn.model_selection import cross_val_score, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import potentials as P
from potentials import HARTREE2CM, BOHR2ANG, v_theta

warnings.filterwarnings("ignore")
T_EVAL = 10.0
CV = KFold(n_splits=5, shuffle=True, random_state=0)
R_GRID_ANG = [2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0]
ANGLES = [(0, 1.0), (90, 0.0), (180, -1.0)]


def logk_at(row, T):
    x = T ** (-1.0 / 6.0)
    return sum(row[f"a{n}"] * x ** n for n in range(6))


def fingerprint_table(lib):
    rows = []
    for _, r in lib.iterrows():
        d = {"id": int(r["id"])}
        for Rang in R_GRID_ANG:
            Rb = Rang / BOHR2ANG
            for adeg, ct in ANGLES:
                v = v_theta(Rb, r["A"], r["a"], r["C6"], r["a1"], r["a2"],
                            r["eta"], ct) * HARTREE2CM
                d[f"V_{int(Rang*10)}_{adeg}"] = np.arcsinh(v / 200.0)
        rows.append(d)
    return pd.DataFrame(rows).set_index("id")


def gbt():
    return HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05,
                                         max_depth=3, random_state=0)


def gp(n):
    k = ConstantKernel() * Matern(length_scale=[1.0] * n, nu=2.5) + WhiteKernel()
    return make_pipeline(StandardScaler(),
                         GaussianProcessRegressor(kernel=k, normalize_y=True,
                                                  n_restarts_optimizer=1,
                                                  random_state=0))


def r2(model, X, y):
    return cross_val_score(model, X, y, cv=CV, scoring="r2").mean()


def main():
    lib = pd.read_csv("potentials_library.csv")
    kt = pd.read_csv("kT_coeffs.csv")

    fp = fingerprint_table(lib)
    fd = P.build_feature_table(weight="depth").add_suffix("_d")
    fn = P.build_feature_table(weight="none").add_suffix("_n")
    base = lib.set_index("id")
    feats = base.join(fd).join(fn).join(fp)
    feats["logB"] = np.log10(feats["B_cm1"])
    feats["log_IV1_d"] = np.log10(np.abs(feats["I_V1_d"]) + 1e-30)
    feats["log_IV1_n"] = np.log10(np.abs(feats["I_V1_n"]) + 1e-30)
    for c in ("eps_0", "eps_90", "eps_180"):
        feats[f"log_{c}"] = np.log10(feats[c])

    FP = list(fp.columns) + ["logB"]            # + mu added per target
    GEOM = ["logB", "log_eps_90", "log_eps_0", "log_eps_180", "d_eps_odd",
            "d_eps_even", "Rm_90", "Rwall_90"]
    DWBA = GEOM + ["log_IV1_d", "I_V1_norm_d", "log_IV1_n", "eta", "a1"]

    def target(partner, ji, jf):
        t = kt[(kt.partner == partner) & (kt.j_i == ji) & (kt.j_f == jf)].copy()
        t["y"] = t.apply(lambda r: logk_at(r, T_EVAL), axis=1)
        t = t[["id", "mu_amu", "y"]].set_index("id")
        d = t.join(feats, how="inner").dropna(subset=["I_V1_d"])
        d["dE"] = (ji * (ji + 1) - jf * (jf + 1)) * d["B_cm1"]
        return d

    # ---- Part A: fingerprint vs derived on He 1->0, both models ----
    d = target("He", 1, 0)
    y = d["y"].values
    print(f"PART A  target = log10 k(He, 1->0, {T_EVAL:.0f}K)   "
          f"N={len(d)}  range[{y.min():.2f},{y.max():.2f}]\n")
    print(f"{'feature set':16} {'#':>3} {'GBT R2':>8} {'GP R2':>8}")
    for name, cols in [("fingerprint", FP + ["mu_amu"]),
                       ("geometry", GEOM + ["mu_amu", "dE"]),
                       ("+ DWBA", DWBA + ["mu_amu", "dE"])]:
        X = d[cols].values
        print(f"{name:16} {len(cols):>3} {r2(gbt(), X, y):8.3f} {r2(gp(len(cols)), X, y):8.3f}")

    # ---- Part B: dj=1 ladder + H2, GBT only ----
    print(f"\nPART B  dj=1 ladder + H2  (GBT R2, target log10 k @ {T_EVAL:.0f}K)\n")
    print(f"{'target':12} {'N':>4} {'fingerprint':>12} {'geometry':>10} {'+DWBA':>8} {'dDWBA':>7}")
    for partner, ji, jf in [("He", 1, 0), ("He", 2, 1), ("He", 3, 2), ("H2", 1, 0)]:
        d = target(partner, ji, jf)
        y = d["y"].values
        rfp = r2(gbt(), d[FP + ["mu_amu"]].values, y)
        rge = r2(gbt(), d[GEOM + ["mu_amu", "dE"]].values, y)
        rdw = r2(gbt(), d[DWBA + ["mu_amu", "dE"]].values, y)
        print(f"{partner} {ji}->{jf:<7} {len(d):>4} {rfp:12.3f} {rge:10.3f} "
              f"{rdw:8.3f} {rdw-rge:+7.3f}")


if __name__ == "__main__":
    main()
