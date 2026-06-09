#!/usr/bin/env python3
"""GP surrogate feasibility test: potential fingerprint -> k(T) profile.

Features (per potential):
  arcsinh(V(R,theta)/S)  for R in {2,3,4,5,6} Ang, theta in {0,90,180} deg  (S=200)
  + log10(B) + mu               -> 17 inputs
Target:
  log10 k(1->0, He, T) at reference temperatures, reconstructed from the fitted
  T^-1/6 coefficients in s2s_kT_coeffs.csv.
5-fold CV; reports R^2 / error per T, a parity plot, and example profiles.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from sklearn.gaussian_process import GaussianProcessRegressor          # noqa: E402
from sklearn.gaussian_process.kernels import (ConstantKernel as C,     # noqa: E402
                                              RBF, WhiteKernel)
from sklearn.preprocessing import StandardScaler                       # noqa: E402
from sklearn.model_selection import KFold                              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
import potential_generator as pg                       # noqa: E402

H = pg.HARTREE2CM
B2A = pg.BOHR2ANG
R_ANG = [2.0, 3.0, 3.5, 4.0, 5.0]       # inner (wall+well); long-range via C6 feature
THETA = {"0": 1.0, "90": 0.0, "180": -1.0}             # cos(theta)
S = 200.0                                              # arcsinh scale [cm^-1]
PARTNER, JI, JF = "He", 1, 0
T_REF = np.array([2.0, 5.0, 10.0, 20.0, 35.0, 50.0])
FIT_DEG = 5


def features(row):
    f = []
    for R in R_ANG:
        Rb = R / B2A
        for ct in THETA.values():
            V = pg.v_theta(Rb, row.A, row.a, row.C6, row.a1, row.a2, row.eta, ct) * H
            f.append(np.arcsinh(V / S))
    return f


def logk_from_coeffs(crow, T):
    x = T ** (-1.0 / 6.0)
    return sum(crow[f"c{n}"] * x ** n for n in range(FIT_DEG + 1))


def main():
    lib = pd.read_csv(os.path.join(REPO, "potentials_library.csv"))
    co = pd.read_csv(os.path.join(HERE, "s2s_kT_coeffs.csv"))
    co = co[(co.partner == PARTNER) & (co.j_i == JI) & (co.j_f == JF)]
    co = co.rename(columns={f"a{n}": f"c{n}" for n in range(FIT_DEG + 1)})
    df = lib.merge(co, on="id", how="inner")
    print(f"training potentials (1->0 {PARTNER}): {len(df)}")

    fcols = [f"V_R{R}_{a}" for R in R_ANG for a in THETA]
    X = np.array([features(r) for _, r in df.iterrows()])
    X = np.column_stack([X, np.log10(df["B_cm1"].values),
                         df["mu_He"].values])
    Y = np.array([[logk_from_coeffs(r, T) for T in T_REF]
                  for _, r in df.iterrows()])
    print(f"features: {X.shape[1]} ({len(fcols)} V(R,theta) + log10B + mu)")

    kernel = (C(1.0) * RBF(length_scale=np.ones(X.shape[1])) +
              WhiteKernel(noise_level=0.05))
    oof = np.zeros_like(Y)
    for tr, te in KFold(5, shuffle=True, random_state=0).split(X):
        sx = StandardScaler().fit(X[tr])
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True,
                                      n_restarts_optimizer=1, alpha=1e-6)
        gp.fit(sx.transform(X[tr]), Y[tr])
        oof[te] = gp.predict(sx.transform(X[te]))

    print("\n  cross-validated accuracy of log10 k (out-of-fold):")
    print(f"  {'T(K)':>6} {'R2':>7} {'MAE(dex)':>9} {'med %err in k':>14}")
    for j, T in enumerate(T_REF):
        yt, yp = Y[:, j], oof[:, j]
        r2 = 1 - np.sum((yt - yp) ** 2) / np.sum((yt - yt.mean()) ** 2)
        mae = np.median(np.abs(yt - yp))
        pct = np.median(10 ** np.abs(yt - yp) - 1) * 100
        print(f"  {T:>6.0f} {r2:>7.3f} {mae:>9.3f} {pct:>13.1f}%")

    # parity at 10 K + example profiles
    j10 = int(np.argmin(np.abs(T_REF - 10)))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    yt, yp = Y[:, j10], oof[:, j10]
    ax[0].scatter(10 ** yt, 10 ** yp, s=9, alpha=0.4)
    lo, hi = 10 ** min(yt.min(), yp.min()), 10 ** max(yt.max(), yp.max())
    ax[0].plot([lo, hi], [lo, hi], "k--", lw=1)
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("true k(1$\\to$0) at 10 K (cm$^3$/s)")
    ax[0].set_ylabel("GP predicted (cm$^3$/s)")
    r2 = 1 - np.sum((yt - yp) ** 2) / np.sum((yt - yt.mean()) ** 2)
    ax[0].set_title(f"parity at 10 K  (R$^2$={r2:.3f}, out-of-fold)")
    ax[0].grid(alpha=.3, which="both")

    rng = np.random.RandomState(1)
    for idx in rng.choice(len(df), 4, replace=False):
        ax[1].plot(T_REF, 10 ** Y[idx], "o-", ms=4, lw=1.2,
                   label=f"pot{int(df.iloc[idx].id)} true")
        ax[1].plot(T_REF, 10 ** oof[idx], "x--", ms=6,
                   color=ax[1].lines[-1].get_color())
    ax[1].set_xlabel("T (K)"); ax[1].set_ylabel("k(1$\\to$0) (cm$^3$/s)")
    ax[1].set_title("example profiles: true (o-) vs GP (x--)")
    ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "gp_demo_1to0_He.png"), dpi=95)
    fig.savefig(os.path.join(HERE, "gp_demo_1to0_He.pdf"))
    print(f"\nwrote gp_demo_1to0_He.png / .pdf")


if __name__ == "__main__":
    main()
