#!/usr/bin/env python3
"""Where does the model fail, and where is the parameter box thin?

Out-of-fold GP residuals (He 1->0, 5 K), then characterise the heavy-tail cases
vs the bulk across physical variables, rank what drives |residual|, and map the
error onto the (B, well-depth) plane to see if hard cases sit in sparse regions.
Outputs hard_region.png + a printed summary to inform the next sampling run.
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
T_EVAL = 5.0
GEOM = ["log_eps_90", "log_eps_0", "log_eps_180", "d_eps_odd", "d_eps_even",
        "Rm_90", "Rwall_90", "logB", "mu_amu", "dE"]


def logk(a, T):
    x = T ** (-1.0 / 6.0)
    return sum(a[n] * x ** n for n in range(6))


def main():
    lib = pd.read_csv("potentials_library.csv").set_index("id")
    kt = pd.read_csv("kT_coeffs.csv")
    t = kt[(kt.partner == "He") & (kt.j_i == 1) & (kt.j_f == 0)].set_index("id")
    g = pd.DataFrame(index=lib.index)
    for c in ("eps_0", "eps_90", "eps_180"):
        g[f"log_{c}"] = np.log10(lib[c])
    for c in ("d_eps_odd", "d_eps_even", "Rm_90", "Rwall_90"):
        g[c] = lib[c]
    g["logB"] = np.log10(lib["B_cm1"])
    d = t.join(g, how="inner")
    d["mu_amu"] = t["mu_amu"]
    d["B"] = lib.loc[d.index, "B_cm1"]; d["dE"] = 2.0 * d["B"]
    d["eps90"] = lib.loc[d.index, "eps_90"]
    d["k_peak"] = t["k_peak"]
    d["y"] = [logk(d.loc[i, [f"a{n}" for n in range(6)]].values, T_EVAL) for i in d.index]

    X = d[GEOM].values; y = d["y"].values
    mu = np.empty(len(y)); sd = np.empty(len(y))
    for tr, te in KFold(5, shuffle=True, random_state=0).split(X):
        sc = StandardScaler().fit(X[tr])
        k = ConstantKernel() * Matern([1.0] * X.shape[1], nu=2.5) + WhiteKernel()
        gp = GaussianProcessRegressor(kernel=k, normalize_y=True, random_state=0)
        gp.fit(sc.transform(X[tr]), y[tr])
        m, s = gp.predict(sc.transform(X[te]), return_std=True)
        mu[te], sd[te] = m, s
    d["resid"] = mu - y
    d["ares"] = np.abs(d["resid"])
    d["depth_thr"] = d["eps90"] / d["dE"]          # well depth / transition gap

    thr = d["ares"].quantile(0.95)
    tail = d[d["ares"] >= thr]; bulk = d[d["ares"] < thr]

    print(f"5 K oof: resid std={d['resid'].std():.3f} dex, "
          f"95th-pct |resid|={thr:.3f} dex, tail N={len(tail)}\n")
    print(f"{'variable':14}{'bulk median':>13}{'tail median':>13}{'ratio':>8}")
    for col, lab in [("B", "B [cm-1]"), ("eps90", "eps90 [cm-1]"),
                     ("d_eps_odd", "d_eps_odd"), ("depth_thr", "eps90/2B"),
                     ("k_peak", "k_peak"), ("mu_amu", "mu")]:
        bm, tm = bulk[col].median(), tail[col].median()
        print(f"{lab:14}{bm:13.4g}{tm:13.4g}{tm/bm:8.2f}")

    imp = HistGradientBoostingRegressor(max_iter=300, max_depth=3, random_state=0)
    from sklearn.inspection import permutation_importance
    imp.fit(X, d["ares"].values)
    pi = permutation_importance(imp, X, d["ares"].values, n_repeats=8, random_state=0)
    order = np.argsort(pi.importances_mean)[::-1]
    print("\ndrivers of |residual| (permutation importance):")
    for i in order[:5]:
        print(f"   {GEOM[i]:12} {pi.importances_mean[i]:.4f}")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    sc0 = ax[0].scatter(d["B"], d["eps90"], c=d["ares"], s=16, cmap="inferno_r",
                        norm=plt.matplotlib.colors.LogNorm(vmin=0.02, vmax=1.0))
    ax[0].scatter(tail["B"], tail["eps90"], s=60, facecolors="none",
                  edgecolors="lime", linewidths=1.3, label="heavy-tail (top 5%)")
    ax[0].set(xscale="log", yscale="log", xlabel="B [cm-1]",
              ylabel="eps_90 [cm-1]", title="error on the (B, depth) plane")
    fig.colorbar(sc0, ax=ax[0], label="|residual| (dex)")
    ax[0].legend(fontsize=8)

    ax[1].scatter(d["depth_thr"], d["ares"], s=14, alpha=0.5, color="purple")
    ax[1].axhline(thr, ls="--", color="grey", label="95th pct")
    ax[1].set(xscale="log", xlabel="eps_90 / 2B  (depth-to-gap)",
              ylabel="|residual| (dex)", title="error vs depth-to-threshold")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig("hard_region.png", dpi=140)
    print("\nsaved hard_region.png")

    # coverage thinness: where are accepted potentials sparse in B?
    print("\nB coverage (accepted potentials per decade):")
    for lo, hi in [(1e-3, 1e-2), (1e-2, 1e-1), (1e-1, 1e0), (1e0, 1e1)]:
        n = ((lib["B_cm1"] >= lo) & (lib["B_cm1"] < hi)).sum()
        print(f"   B in [{lo:.0e},{hi:.0e}): {n:4d}")


if __name__ == "__main__":
    main()
