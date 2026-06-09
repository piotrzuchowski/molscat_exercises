#!/usr/bin/env python3
"""Out-of-fold validation figure for the GP (geometry + DWBA) on He 1->0.

Manual 5-fold CV collecting held-out predictions AND GP sigma, then:
  (a) predicted vs actual log10 k, coloured by predicted sigma
  (b) residual histogram (dex)
  (c) calibration: z = residual/sigma vs standard normal, coverage
  (d) |residual| vs predicted sigma (does the GP know when it's unsure?)
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

import potentials as P

warnings.filterwarnings("ignore")
T_EVAL = 5.0


def logk_at(row, T):
    x = T ** (-1.0 / 6.0)
    return sum(row[f"a{n}"] * x ** n for n in range(6))


def main():
    lib = pd.read_csv("potentials_library.csv")
    kt = pd.read_csv("kT_coeffs.csv")
    fd = P.build_feature_table(weight="depth").add_suffix("_d")
    fn = P.build_feature_table(weight="none").add_suffix("_n")
    feats = lib.set_index("id").join(fd).join(fn)
    feats["logB"] = np.log10(feats["B_cm1"])
    feats["log_IV1_d"] = np.log10(np.abs(feats["I_V1_d"]) + 1e-30)
    feats["log_IV1_n"] = np.log10(np.abs(feats["I_V1_n"]) + 1e-30)
    for c in ("eps_0", "eps_90", "eps_180"):
        feats[f"log_{c}"] = np.log10(feats[c])

    t = kt[(kt.partner == "He") & (kt.j_i == 1) & (kt.j_f == 0)].copy()
    t["y"] = t.apply(lambda r: logk_at(r, T_EVAL), axis=1)
    t = t[["id", "mu_amu", "y"]].set_index("id")
    d = t.join(feats, how="inner").dropna(subset=["I_V1_d"])
    d["dE"] = 2.0 * d["B_cm1"]

    cols = ["logB", "log_eps_90", "log_eps_0", "log_eps_180", "d_eps_odd",
            "d_eps_even", "Rm_90", "Rwall_90", "log_IV1_d", "I_V1_norm_d",
            "log_IV1_n", "eta", "a1", "mu_amu", "dE"]
    X = d[cols].values
    y = d["y"].values
    n = len(y)

    mu = np.empty(n); sd = np.empty(n); fold_r2 = []
    cv = KFold(n_splits=5, shuffle=True, random_state=0)
    for tr, te in cv.split(X):
        sc = StandardScaler().fit(X[tr])
        k = ConstantKernel() * Matern([1.0] * X.shape[1], nu=2.5) + WhiteKernel()
        gp = GaussianProcessRegressor(kernel=k, normalize_y=True,
                                      n_restarts_optimizer=1, random_state=0)
        gp.fit(sc.transform(X[tr]), y[tr])
        m, s = gp.predict(sc.transform(X[te]), return_std=True)
        mu[te], sd[te] = m, s
        fold_r2.append(1 - np.sum((m - y[te]) ** 2) / np.sum((y[te] - y[te].mean()) ** 2))

    res = mu - y
    rmse = np.sqrt(np.mean(res ** 2))
    r2 = 1 - np.sum(res ** 2) / np.sum((y - y.mean()) ** 2)
    z = res / sd
    cov1 = np.mean(np.abs(z) < 1.0)
    cov2 = np.mean(np.abs(z) < 2.0)

    fig, ax = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Out-of-fold validation: GP (geometry+DWBA), "
                 f"log10 k(He, 1->0, {T_EVAL:.0f} K), N={n}", fontsize=13)

    # (a) predicted vs actual
    a = ax[0, 0]
    sc = a.scatter(y, mu, c=sd, cmap="viridis", s=14, alpha=0.8)
    lim = [min(y.min(), mu.min()) - 0.1, max(y.max(), mu.max()) + 0.1]
    a.plot(lim, lim, "k--", lw=1)
    a.set(xlim=lim, ylim=lim, xlabel="actual log10 k", ylabel="predicted log10 k")
    a.set_title(f"predicted vs actual   R²={r2:.3f}, RMSE={rmse:.3f} dex")
    fig.colorbar(sc, ax=a, label="GP σ (dex)", fraction=0.046)

    # (b) residual histogram
    b = ax[0, 1]
    b.hist(res, bins=40, color="steelblue", edgecolor="k", alpha=0.8)
    b.axvline(0, color="k", lw=1)
    b.set(xlabel="residual = pred - actual (dex)", ylabel="count")
    b.set_title(f"residuals   mean={res.mean():+.3f}, std={res.std():.3f} dex, "
                f"95%|{np.percentile(np.abs(res),95):.2f}")

    # (c) calibration: z-score vs N(0,1)
    c = ax[1, 0]
    c.hist(z, bins=40, density=True, color="salmon", edgecolor="k", alpha=0.8)
    xx = np.linspace(-4, 4, 200)
    c.plot(xx, norm.pdf(xx), "k-", lw=2, label="N(0,1) ideal")
    c.set(xlabel="z = residual / GP σ", ylabel="density")
    c.set_title(f"calibration   coverage ±1σ={cov1:.2f} (0.68), ±2σ={cov2:.2f} (0.95)")
    c.legend(fontsize=8)

    # (d) |residual| vs predicted sigma
    dax = ax[1, 1]
    dax.scatter(sd, np.abs(res), s=14, alpha=0.6, color="purple")
    m = max(sd.max(), np.abs(res).max())
    dax.plot([0, m], [0, m], "k--", lw=1, label="|res| = σ")
    dax.set(xlabel="GP σ (dex)", ylabel="|residual| (dex)")
    dax.set_title("does the GP know when it's unsure?")
    dax.legend(fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fname = f"validation_He{int(T_EVAL)}.png"
    fig.savefig(fname, dpi=140)
    print(f"oof R²={r2:.3f}  RMSE={rmse:.3f} dex  fold R²={np.round(fold_r2,3)}")
    print(f"residual std={res.std():.3f} dex  95%={np.percentile(np.abs(res),95):.2f} dex")
    print(f"coverage ±1σ={cov1:.2f} (ideal 0.68)  ±2σ={cov2:.2f} (ideal 0.95)")
    print(f"saved {fname}")


if __name__ == "__main__":
    main()
