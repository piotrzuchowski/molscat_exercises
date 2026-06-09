#!/usr/bin/env python3
"""Learning curve: is the ~0.9 ceiling data-limited or information-limited?

Fixed held-out test set; train GBT and GP on growing random subsets of the rest;
R^2 on the fixed test vs train size. Flat at the top => more data won't help
(information/physics ceiling). Still rising => more data helps.
Target: log10 k(He, 1->0, 10 K), geometry features.
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
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

warnings.filterwarnings("ignore")
T_EVAL = 10.0
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
    d["mu_amu"] = t["mu_amu"]; d["dE"] = 2.0 * lib.loc[d.index, "B_cm1"]
    d["y"] = [logk(d.loc[i, [f"a{n}" for n in range(6)]].values, T_EVAL) for i in d.index]

    rng = np.random.default_rng(0)
    idx = rng.permutation(len(d))
    test = idx[:200]; trainpool = idx[200:]
    X = d[GEOM].values; y = d["y"].values
    Xte, yte = X[test], y[test]

    sizes = [50, 100, 175, 300, 450, 600, len(trainpool)]
    rows = []
    for nsz in sizes:
        rg, rp = [], []
        for rep in range(6):
            sub = rng.choice(trainpool, size=nsz, replace=False)
            gbt = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05,
                                                max_depth=3, random_state=rep)
            gbt.fit(X[sub], y[sub])
            rg.append(r2_score(yte, gbt.predict(Xte)))
            if nsz <= 600:                      # GP only where cheap
                sc = StandardScaler().fit(X[sub])
                k = ConstantKernel() * Matern([1.0] * X.shape[1], nu=2.5) + WhiteKernel()
                gp = GaussianProcessRegressor(kernel=k, normalize_y=True, random_state=rep)
                gp.fit(sc.transform(X[sub]), y[sub])
                rp.append(r2_score(yte, gp.predict(sc.transform(Xte))))
        rows.append({"n": nsz, "gbt_mean": np.mean(rg), "gbt_std": np.std(rg),
                     "gp_mean": np.mean(rp) if rp else np.nan,
                     "gp_std": np.std(rp) if rp else 0.0})
        print(f"  N={nsz:4d}  GBT R2={np.mean(rg):.3f}+-{np.std(rg):.3f}"
              + (f"  GP R2={np.mean(rp):.3f}" if rp else ""))

    r = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(r.n, r.gbt_mean, yerr=r.gbt_std, marker="o", capsize=3,
                label="GBT", color="C0")
    m = r.gp_mean.notna()
    ax.errorbar(r.n[m], r.gp_mean[m], yerr=r.gp_std[m], marker="s", capsize=3,
                label="GP", color="C1")
    ax.axhline(0.90, ls="--", color="grey", lw=1, label="current ~0.90")
    ax.set_xscale("log")
    ax.set(xlabel="training potentials", ylabel="R$^2$ on fixed 200-pt test",
           title="Learning curve: log10 k(He, 1->0, 10 K)")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig("learning_curve.png", dpi=140)
    print("saved learning_curve.png")
    # slope diagnostic over last decade
    last = r.iloc[-1]; prev = r.iloc[-3]
    print(f"\nslope (last leg): dR2={last.gbt_mean-prev.gbt_mean:+.3f} over "
          f"N {int(prev.n)}->{int(last.n)} ({last.n/prev.n:.1f}x data)")


if __name__ == "__main__":
    main()
