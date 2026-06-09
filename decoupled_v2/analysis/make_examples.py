#!/usr/bin/env python3
"""Example k(T) curves: true (dataset fit) vs GP prediction with uncertainty.

10 potentials spanning large/small V0 (well depth), V1 (odd anisotropy d_eps_odd)
and V2 (even anisotropy d_eps_even), shown for transitions 1->0, 3->2 (dj=1, V1)
and 2->0 (dj=2, V2), He partner.

The GP (geometry features) is trained per temperature node on ALL potentials
EXCEPT the 10 examples, so every plotted prediction is genuinely out-of-sample.
Band = +-1 GP sigma in log10 k.
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

GEOM = ["log_eps_90", "log_eps_0", "log_eps_180", "d_eps_odd", "d_eps_even",
        "Rm_90", "Rwall_90", "logB"]
COLS = GEOM + ["mu_amu", "dE"]
TRANS = [("He", 1, 0, "1->0"), ("He", 3, 2, "3->2"), ("He", 2, 0, "2->0")]
T_NODE = np.geomspace(5.0, 50.0, 9)
T_FINE = np.linspace(5.0, 50.0, 60)


def logk(a, T):
    x = T ** (-1.0 / 6.0)
    return sum(a[n] * x ** n for n in range(6))


def farthest_points(P, k, seed=0):
    """Farthest-point sampling: k diverse rows of standardized array P."""
    P = (P - P.mean(0)) / (P.std(0) + 1e-9)
    idx = [int(np.argmax(np.linalg.norm(P - P.mean(0), axis=1)))]
    for _ in range(k - 1):
        d = np.min([np.linalg.norm(P - P[i], axis=1) for i in idx], axis=0)
        idx.append(int(np.argmax(d)))
    return idx


def main():
    lib = pd.read_csv("potentials_library.csv").set_index("id")
    kt = pd.read_csv("kT_coeffs.csv")
    g = pd.DataFrame(index=lib.index)
    for c in ("eps_0", "eps_90", "eps_180"):
        g[f"log_{c}"] = np.log10(lib[c])
    for c in ("d_eps_odd", "d_eps_even", "Rm_90", "Rwall_90", "B_cm1"):
        g[c] = lib[c]
    g["logB"] = np.log10(lib["B_cm1"])

    def tdf(partner, ji, jf):
        t = kt[(kt.partner == partner) & (kt.j_i == ji) & (kt.j_f == jf)].copy()
        t = t.set_index("id")
        d = t.join(g, how="inner")
        d["dE"] = (ji * (ji + 1) - jf * (jf + 1)) * d["B_cm1"]
        return d

    dfs = {lab: tdf(p, i, f) for p, i, f, lab in TRANS}
    common = set.intersection(*[set(d.index) for d in dfs.values()])
    common = sorted(common)
    sel_feats = lib.loc[common, ["eps_90", "d_eps_odd", "d_eps_even"]].copy()
    sel_feats["eps_90"] = np.log10(sel_feats["eps_90"])
    sel_feats["d_eps_odd"] = np.abs(sel_feats["d_eps_odd"])
    sel_feats["d_eps_even"] = np.abs(sel_feats["d_eps_even"])
    pick = [common[i] for i in farthest_points(sel_feats.values, 10)]
    print("selected ids:", pick)

    # train per-transition, per-node GP on all-but-examples; predict examples
    pred = {lab: {cid: {"m": [], "s": []} for cid in pick} for _, _, _, lab in TRANS}
    for (_, _, _, lab) in TRANS:
        d = dfs[lab]
        tr = d.drop(index=pick)
        Xtr = tr[COLS].values
        sc = StandardScaler().fit(Xtr)
        Xex = sc.transform(d.loc[pick, COLS].values)
        for Tn in T_NODE:
            ytr = np.array([logk(tr.loc[i, [f"a{n}" for n in range(6)]].values, Tn)
                            for i in tr.index])
            k = ConstantKernel() * Matern([1.0] * len(COLS), nu=2.5) + WhiteKernel()
            gp = GaussianProcessRegressor(kernel=k, normalize_y=True,
                                          n_restarts_optimizer=1, random_state=0)
            gp.fit(sc.transform(Xtr), ytr)
            m, s = gp.predict(Xex, return_std=True)
            for j, cid in enumerate(pick):
                pred[lab][cid]["m"].append(m[j])
                pred[lab][cid]["s"].append(s[j])
        print(f"  trained {lab}")

    fig, axes = plt.subplots(10, 3, figsize=(13, 26), sharex=True)
    fig.suptitle("k(T): dataset fit (black) vs GP prediction +- 1 sigma (blue), "
                 "out-of-sample He", fontsize=14)
    for r, cid in enumerate(pick):
        e90 = lib.loc[cid, "eps_90"]; dodd = lib.loc[cid, "d_eps_odd"]
        deven = lib.loc[cid, "d_eps_even"]
        for c, (_, _, _, lab) in enumerate(TRANS):
            ax = axes[r, c]
            a = dfs[lab].loc[cid, [f"a{n}" for n in range(6)]].values
            ax.plot(T_FINE, 10 ** logk(a, T_FINE), "k-", lw=1.6, label="dataset")
            m = np.array(pred[lab][cid]["m"]); s = np.array(pred[lab][cid]["s"])
            ax.plot(T_NODE, 10 ** m, "b--o", ms=3, lw=1.2, label="GP mean")
            ax.fill_between(T_NODE, 10 ** (m - s), 10 ** (m + s),
                            color="b", alpha=0.2, label="GP +-1 sigma")
            ax.set_yscale("log")
            if r == 0:
                ax.set_title(lab, fontsize=12)
            if c == 0:
                ax.set_ylabel(f"id {cid}\n$\\epsilon_{{90}}$={e90:.0f}\n"
                              f"$\\Delta\\epsilon_{{odd}}$={dodd:.0f}\n"
                              f"$\\Delta\\epsilon_{{even}}$={deven:.0f}",
                              fontsize=8)
            if r == 9:
                ax.set_xlabel("T [K]")
            ax.tick_params(labelsize=7)
            if r == 0 and c == 2:
                ax.legend(fontsize=7, loc="best")
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig("example_kT_curves.png", dpi=120)
    print("saved example_kT_curves.png")


if __name__ == "__main__":
    main()
