#!/usr/bin/env python3
"""Surrogate with decorrelated (orthonormal-basis) targets.

features -> c_m (Chebyshev-orthonormal coefficients of ln k over a global T grid)
-> rebuild k(T). Compares against predicting on the raw grid, and reports the
result back in the a_n monomial convention. CV over the 1->0 He curves.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (ConstantKernel as C, RBF,
                                              Matern, WhiteKernel)
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold

import gp_demo as g                       # reuse feature builder
from logk_basis import LogKBasis

HERE = g.HERE
REPO = g.REPO
T = np.linspace(5.0, 50.0, 15)
basis = LogKBasis(T, deg=5)
LN10 = np.log(10.0)


def main():
    lib = pd.read_csv(os.path.join(REPO, "potentials_library.csv"))
    co = pd.read_csv(os.path.join(HERE, "s2s_kT_coeffs.csv"))
    co = co[(co.partner == "He") & (co.j_i == 1) & (co.j_f == 0)]
    co = co.rename(columns={f"a{n}": f"c{n}" for n in range(6)})
    df = lib.merge(co, on="id", how="inner").reset_index(drop=True)
    print(f"curves (1->0 He): {len(df)}; global grid {T[0]:.0f}-{T[-1]:.0f} K x{len(T)}")

    # features (17) and ln k curves on the global grid (from saved log10 fit)
    X = np.array([g.features(r) for _, r in df.iterrows()])
    X = np.column_stack([X, np.log10(df["C6"].values),       # long-range tail
                         np.log10(df["B_cm1"].values), df["mu_He"].values])
    x = T ** (-1.0 / 6.0)
    lnk = np.array([[(sum(r[f"c{n}"] * xx ** n for n in range(6))) * LN10
                     for xx in x] for _, r in df.iterrows()])
    Cmat = basis.fit_many(lnk)            # (n, 6) decorrelated targets

    kernel = (C(1.0) * Matern(np.ones(X.shape[1]), nu=2.5) + WhiteKernel(0.05))
    oofC = np.zeros_like(Cmat)
    for tr, te in KFold(5, shuffle=True, random_state=0).split(X):
        sx = StandardScaler().fit(X[tr])
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True,
                                      n_restarts_optimizer=0, alpha=1e-6)
        gp.fit(sx.transform(X[tr]), Cmat[tr])     # multi-output over modes
        oofC[te] = gp.predict(sx.transform(X[te]))

    print("\n  per-mode out-of-fold R^2 (variance share):")
    for m in range(6):
        yt, yp = Cmat[:, m], oofC[:, m]
        r2 = 1 - np.sum((yt - yp) ** 2) / np.sum((yt - yt.mean()) ** 2)
        share = Cmat[:, m].var() / Cmat.var(0).sum()
        print(f"    c{m}:  R2 = {r2:6.3f}   ({share*100:5.1f}% of variance)")

    # rebuild curves and score on k
    lnk_pred = basis.eval(oofC)                  # (n, nT)
    dex = np.abs(lnk_pred - lnk) / LN10          # error in log10 k
    pct = np.median(np.expm1(np.abs(lnk_pred - lnk))) * 100
    r2c = 1 - np.sum((lnk_pred - lnk) ** 2) / np.sum((lnk - lnk.mean()) ** 2)
    print(f"\n  reconstructed curve: overall R^2(ln k) = {r2c:.3f}, "
          f"median |err| = {np.median(dex):.3f} dex, median %err in k = {pct:.1f}%")

    # report one prediction back in a_n (natural-log monomial) convention
    a0 = basis.to_monomial(oofC[0])
    print(f"\n  example predicted a_n (ln k = sum a_n T^-n/6) for pot{int(df.id[0])}:")
    print("   ", np.array2string(a0, precision=3, suppress_small=False))

    # plots: parity of c0 (level) + example rebuilt profiles
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    ax[0].scatter(Cmat[:, 0], oofC[:, 0], s=9, alpha=0.4)
    lo, hi = Cmat[:, 0].min(), Cmat[:, 0].max()
    ax[0].plot([lo, hi], [lo, hi], "k--", lw=1)
    r2 = 1 - np.sum((Cmat[:, 0] - oofC[:, 0]) ** 2) / np.sum(
        (Cmat[:, 0] - Cmat[:, 0].mean()) ** 2)
    ax[0].set_xlabel("true c$_0$ (level mode)")
    ax[0].set_ylabel("GP predicted c$_0$")
    ax[0].set_title(f"level mode c$_0$  (R$^2$={r2:.3f}, 98.5% of variance)")
    ax[0].grid(alpha=.3)
    Tf = np.linspace(5, 50, 100)
    rng = np.random.RandomState(2)
    for idx in rng.choice(len(df), 4, replace=False):
        kt = basis.k_of_T(Cmat[idx], Tf)
        kp = basis.k_of_T(oofC[idx], Tf)
        ax[1].plot(Tf, kt, "-", lw=1.3, label=f"pot{int(df.id[idx])}")
        ax[1].plot(Tf, kp, "--", color=ax[1].lines[-1].get_color())
    ax[1].set_xlabel("T (K)"); ax[1].set_ylabel("k(1$\\to$0) (cm$^3$/s)")
    ax[1].set_title("rebuilt curves: true (-) vs GP (--)")
    ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "gp_demo2_chebyshev.png"), dpi=95)
    fig.savefig(os.path.join(HERE, "gp_demo2_chebyshev.pdf"))
    print("\n  wrote gp_demo2_chebyshev.png / .pdf")


if __name__ == "__main__":
    main()
