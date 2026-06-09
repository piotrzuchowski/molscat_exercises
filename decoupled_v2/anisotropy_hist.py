#!/usr/bin/env python3
"""Variability check: V1/V0 and V2/V0 evaluated at the V0 minimum, over the
whole decoupled_v2 library.

  V0(R) =      A e^{-aR} -      C6 R^{-6}
  V1(R) =  eta A e^{-aR} - a1   C6 R^{-6}
  V2(R) = zeta A e^{-aR} - a2   C6 R^{-6}

Ratios are taken at Rm = argmin V0 (the well bottom). Histograms + summary.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
lib = pd.read_csv(os.path.join(HERE, "potentials_library.csv"))


def v0(R, A, a, C6):
    return A * np.exp(-a * R) - C6 * R ** -6


def comp(R, A, a, C6, w_exp, w_disp):     # generic Legendre component
    return w_exp * A * np.exp(-a * R) - w_disp * C6 * R ** -6


r1, r2, rm = [], [], []
for _, p in lib.iterrows():
    res = minimize_scalar(lambda R: v0(R, p.A, p.a, p.C6),
                          bounds=(2.5, 40.0), method="bounded")
    Rm = res.x
    V0 = v0(Rm, p.A, p.a, p.C6)
    V1 = comp(Rm, p.A, p.a, p.C6, p.eta, p.a1)
    V2 = comp(Rm, p.A, p.a, p.C6, p.zeta, p.a2)
    r1.append(V1 / V0); r2.append(V2 / V0); rm.append(Rm * 0.529177210903)

lib["V1_over_V0"] = r1
lib["V2_over_V0"] = r2
lib["Rm_iso_ang"] = rm
lib[["id", "V1_over_V0", "V2_over_V0", "Rm_iso_ang"]].to_csv(
    os.path.join(HERE, "anisotropy_at_min.csv"), index=False)

r1 = np.array(r1); r2 = np.array(r2)
print(f"N = {len(r1)} potentials,  ratios at V0 minimum")
for name, r in (("V1/V0", r1), ("V2/V0", r2)):
    q = np.percentile(r, [1, 25, 50, 75, 99])
    print(f"  {name}: mean={r.mean():+.3f} std={r.std():.3f} "
          f"min={r.min():+.3f} max={r.max():+.3f}  "
          f"[1,25,50,75,99]%={np.array2string(q, precision=2, sign='+')}")
print(f"  corr(V1/V0, V2/V0) = {np.corrcoef(r1, r2)[0,1]:+.3f}")

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
for a, r, name, c in ((ax[0], r1, "V1/V0", "tab:blue"),
                      (ax[1], r2, "V2/V0", "tab:red")):
    a.hist(r, bins=50, color=c, alpha=0.8, edgecolor="k", lw=0.3)
    a.axvline(0, color="k", lw=0.8, ls="--")
    a.set_xlabel(f"{name}  at  R$_m$(V0)")
    a.set_ylabel("count")
    a.set_title(f"{name}   (std={r.std():.2f}, range {r.min():+.2f}..{r.max():+.2f})")
    a.grid(alpha=0.3)
ax[2].scatter(r1, r2, s=8, alpha=0.4)
ax[2].axhline(0, color="k", lw=0.6); ax[2].axvline(0, color="k", lw=0.6)
ax[2].set_xlabel("V1/V0"); ax[2].set_ylabel("V2/V0")
ax[2].set_title(f"joint coverage (corr={np.corrcoef(r1, r2)[0,1]:+.2f})")
ax[2].grid(alpha=0.3)
fig.suptitle("decoupled_v2 anisotropy at the V0 minimum (946 potentials)",
             fontsize=12, weight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(HERE, "anisotropy_hist.png"), dpi=95)
fig.savefig(os.path.join(HERE, "anisotropy_hist.pdf"))
print("wrote anisotropy_hist.png/.pdf and anisotropy_at_min.csv")
