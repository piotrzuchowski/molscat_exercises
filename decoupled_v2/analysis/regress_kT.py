#!/usr/bin/env python3
"""Regress 1->0 He rate coefficient on DWBA overlap features + physics.

Target: log10 k(10 K), reconstructed from the a0..a5  (log10 k = sum a_n T^-n/6).
Feature ablation: kinematics -> + geometry -> + DWBA overlap, with 5-fold CV R^2,
to compare against the README's R^2 ~ 0.5 fingerprint baseline.
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

warnings.filterwarnings("ignore")
T_EVAL = 10.0
CV = KFold(n_splits=5, shuffle=True, random_state=0)


def logk_at(row, T):
    x = T ** (-1.0 / 6.0)
    return sum(row[f"a{n}"] * x ** n for n in range(6))


def main():
    lib = pd.read_csv("potentials_library.csv")
    kt = pd.read_csv("kT_coeffs.csv")

    # DWBA overlap features (two weightings) keyed by id
    f_depth = P.build_feature_table(weight="depth").add_suffix("_d")
    f_none = P.build_feature_table(weight="none").add_suffix("_n")
    feats = lib.set_index("id").join(f_depth).join(f_none)

    tgt = kt[(kt.partner == "He") & (kt.j_i == 1) & (kt.j_f == 0)].copy()
    tgt["y"] = tgt.apply(lambda r: logk_at(r, T_EVAL), axis=1)
    tgt = tgt[["id", "mu_amu", "k_peak", "y"]].set_index("id")  # drop a0..a5 (clash)
    d = tgt.join(feats, how="inner").dropna(subset=["I_V1_d"])
    print(f"{len(d)} potentials with 1->0 He data + features\n")

    # engineered columns
    d["logB"] = np.log10(d["B_cm1"])
    d["dE"] = 2.0 * d["B_cm1"]                       # 1->0 threshold spacing
    d["log_IV1_d"] = np.log10(np.abs(d["I_V1_d"]) + 1e-30)
    d["log_IV1_n"] = np.log10(np.abs(d["I_V1_n"]) + 1e-30)
    d["log_eps90"] = np.log10(d["eps_90"])
    d["log_eps0"] = np.log10(d["eps_0"])

    sets = {
        "kinematics      ": ["logB", "mu_amu", "dE"],
        "+ geometry      ": ["logB", "mu_amu", "dE", "log_eps90", "log_eps0",
                             "d_eps_odd", "Rm_90", "Rwall_90"],
        "+ DWBA overlap  ": ["logB", "mu_amu", "dE", "log_eps90", "log_eps0",
                             "d_eps_odd", "Rm_90", "Rwall_90",
                             "log_IV1_d", "I_V1_norm_d", "log_IV1_n", "eta", "a1"],
        "DWBA-lean       ": ["logB", "mu_amu", "log_eps90", "log_IV1_d",
                             "I_V1_norm_d", "eta", "a1"],
    }

    y = d["y"].values
    print(f"target = log10 k(1->0, He, {T_EVAL:.0f} K)   "
          f"[range {y.min():.2f} .. {y.max():.2f}]\n")
    print(f"{'feature set':18} {'#':>3} {'GBT R2':>8} {'GP R2':>8}")
    for name, cols in sets.items():
        X = d[cols].values
        gbt = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05,
                                            max_depth=3, random_state=0)
        r_gbt = cross_val_score(gbt, X, y, cv=CV, scoring="r2").mean()
        k = (ConstantKernel() * Matern(length_scale=[1.0] * X.shape[1], nu=2.5)
             + WhiteKernel())
        gp = make_pipeline(StandardScaler(),
                           GaussianProcessRegressor(kernel=k, normalize_y=True,
                                                    n_restarts_optimizer=1,
                                                    random_state=0))
        r_gp = cross_val_score(gp, X, y, cv=CV, scoring="r2").mean()
        print(f"{name} {len(cols):>3} {r_gbt:8.3f} {r_gp:8.3f}")

    print("\n(README fingerprint baseline on this target: R2 ~ 0.5)")


if __name__ == "__main__":
    main()
