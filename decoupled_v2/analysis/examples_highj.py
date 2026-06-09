#!/usr/bin/env python3
"""Cleaner k(T) example figure for HIGHER-j transitions (4->3, 6->5, 8->7).

GP trained on 600 potentials (O(N^3): 600 vs 900 ~3x cheaper, ~same accuracy).
6 example potentials, held out of training -> out-of-sample. Black = dataset
fit; blue = GP mean; band = +-1 sigma.
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
N_TRAIN = 600
N_CASES = 6
GEOM = ["log_eps_90", "log_eps_0", "log_eps_180", "d_eps_odd", "d_eps_even",
        "Rm_90", "Rwall_90", "logB", "mu_amu", "dE"]
TRANS = [("He", 5, 4, "5->4"), ("He", 6, 5, "6->5"), ("He", 7, 6, "7->6")]
T_NODE = np.geomspace(5.0, 50.0, 9)
T_FINE = np.linspace(5.0, 50.0, 60)
RNG = np.random.default_rng(1)


def logk(a, T):
    x = T ** (-1.0 / 6.0)
    return sum(a[n] * x ** n for n in range(6))


def fps(P, k):
    P = (P - P.mean(0)) / (P.std(0) + 1e-9)
    idx = [int(np.argmax(np.linalg.norm(P - P.mean(0), axis=1)))]
    for _ in range(k - 1):
        dmin = np.min([np.linalg.norm(P - P[i], axis=1) for i in idx], axis=0)
        idx.append(int(np.argmax(dmin)))
    return idx


def main():
    lib = pd.read_csv("potentials_library.csv").set_index("id")
    kt = pd.read_csv("kT_coeffs.csv")
    g = pd.DataFrame(index=lib.index)
    for c in ("eps_0", "eps_90", "eps_180"):
        g[f"log_{c}"] = np.log10(lib[c])
    for c in ("d_eps_odd", "d_eps_even", "Rm_90", "Rwall_90"):
        g[c] = lib[c]
    g["logB"] = np.log10(lib["B_cm1"])

    def tdf(p, i, f):
        t = kt[(kt.partner == p) & (kt.j_i == i) & (kt.j_f == f)].set_index("id")
        d = t.join(g, how="inner")
        d["mu_amu"] = t["mu_amu"]
        d["dE"] = (i * (i + 1) - f * (f + 1)) * lib.loc[d.index, "B_cm1"]
        return d

    dfs = {lab: tdf(p, i, f) for p, i, f, lab in TRANS}
    common = sorted(set.intersection(*[set(d.index) for d in dfs.values()]))
    sf = lib.loc[common, ["eps_90", "d_eps_odd", "B_cm1"]].copy()
    sf["eps_90"] = np.log10(sf["eps_90"]); sf["B_cm1"] = np.log10(sf["B_cm1"])
    sf["d_eps_odd"] = np.abs(sf["d_eps_odd"])
    pick = [common[i] for i in fps(sf.values, N_CASES)]
    print("cases:", pick, f"({len(common)} common potentials)")

    pred = {lab: {c: {"m": [], "s": []} for c in pick} for _, _, _, lab in TRANS}
    for (_, _, _, lab) in TRANS:
        d = dfs[lab]
        pool = [i for i in d.index if i not in pick]
        tr = d.loc[RNG.choice(pool, size=min(N_TRAIN, len(pool)), replace=False)]
        Xtr = tr[GEOM].values
        sc = StandardScaler().fit(Xtr)
        Xex = sc.transform(d.loc[pick, GEOM].values)
        for Tn in T_NODE:
            ytr = np.array([logk(tr.loc[i, [f"a{n}" for n in range(6)]].values, Tn)
                            for i in tr.index])
            k = ConstantKernel() * Matern([1.0] * len(GEOM), nu=2.5) + WhiteKernel()
            gp = GaussianProcessRegressor(kernel=k, normalize_y=True, random_state=0)
            gp.fit(sc.transform(Xtr), ytr)
            m, s = gp.predict(Xex, return_std=True)
            for j, c in enumerate(pick):
                pred[lab][c]["m"].append(m[j]); pred[lab][c]["s"].append(s[j])
        print(f"  trained {lab} on {len(tr)}")

    fig, ax = plt.subplots(N_CASES, 3, figsize=(12, 2.4 * N_CASES), sharex=True)
    fig.suptitle("Higher-j k(T): dataset (black) vs GP +-1 sigma (blue), "
                 "out-of-sample He  [GP trained on 600]", fontsize=13)
    for r, cid in enumerate(pick):
        e90 = lib.loc[cid, "eps_90"]; dodd = lib.loc[cid, "d_eps_odd"]
        B = lib.loc[cid, "B_cm1"]
        for c, (_, _, _, lab) in enumerate(TRANS):
            a = ax[r, c]
            av = dfs[lab].loc[cid, [f"a{n}" for n in range(6)]].values
            a.plot(T_FINE, 10 ** logk(av, T_FINE), "k-", lw=1.8, label="dataset")
            m = np.array(pred[lab][cid]["m"]); s = np.array(pred[lab][cid]["s"])
            a.plot(T_NODE, 10 ** m, "-", color="C0", lw=1.6, label="GP mean")
            a.fill_between(T_NODE, 10 ** (m - s), 10 ** (m + s), color="C0",
                           alpha=0.25, label="GP $\\pm1\\sigma$")
            a.set_yscale("log")
            if r == 0:
                a.set_title(lab, fontsize=12)
            if c == 0:
                a.set_ylabel(f"id {cid}\nB={B:.2g}\n$\\epsilon_{{90}}$={e90:.0f}\n"
                             f"$\\Delta\\epsilon_{{odd}}$={dodd:.0f}", fontsize=8)
            if r == N_CASES - 1:
                a.set_xlabel("T [K]")
            a.tick_params(labelsize=8)
            if r == 0 and c == 2:
                a.legend(fontsize=8, loc="best")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig("examples_highj.png", dpi=135)
    print("saved examples_highj.png")


if __name__ == "__main__":
    main()
