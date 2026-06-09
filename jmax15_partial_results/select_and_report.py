#!/usr/bin/env python3
"""Select sample potentials by (V0 depth, V1/V2 anisotropy) and make a k(T) PDF.

Cases (anisotropy x depth), B in 0.1-1 cm^-1, from COMPLETED potentials only:
  1 = V1 & V2 weak     2 = V1 large, V2 weak     3 = V1 weak, V2 large
  A = V0 ~ 20 cm^-1    B = V0 ~ 40 cm^-1
-> 1A 1B 2A 2B 3A 3B.

Definitions:
  V0 depth  = isotropic well depth of V0(R)=A e^{-aR}-C6/R^6  (eps0_iso)
  V1 strength = |d_eps_odd| / eps0_iso    (odd, lambda=1 contrast)
  V2 strength = |d_eps_even| / eps0_iso   (even, lambda=2 contrast)
Best-match selection (closest available; actual values reported).

For each case: de-excitation k(i->f, T) for He and H2, markers = directly
computed points (k_of_T at 25 T), continuous line = log10 k = sum a_n (T^-1/6)^n
fit. Linear scales. Multi-page PDF.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt              # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "smoke_test"))
import potential_generator as pg             # noqa: E402
import energy_kinteg_helper as ek            # noqa: E402
from parse_xs import parse_file              # noqa: E402

PROD = os.path.join(REPO, "production")
H = pg.HARTREE2CM
T_MARK = np.geomspace(2.0, 50.0, 25)         # calculated points
T_FINE = np.linspace(2.0, 50.0, 200)         # smooth fit line
FIT_DEG = 5
KMIN = 1.0e-12
# fixed low-j set: Dj=-1 ladder (V1/odd-driven) + Dj=-2 (V2/even-driven)
FIXED_TRANS = [(1, 0), (2, 1), (3, 2), (2, 0), (3, 1)]
PDF = os.path.join(HERE, "kT_cases_report.pdf")

# (case, mode, depth target)
CASES = [("1A", "weakweak", 20), ("1B", "weakweak", 40),
         ("2A", "v1large", 20), ("2B", "v1large", 40),
         ("3A", "v2large", 20), ("3B", "v2large", 40)]


def complete_ids(man):
    def done(p):
        try:
            return "This run used" in open(os.path.join(PROD, p)).read()
        except OSError:
            return False
    m = man.copy()
    m["out"] = m.apply(lambda r: os.path.join(
        r["dir"], r["input"].replace(".input", ".out")), axis=1)
    m["done"] = m["out"].apply(done)
    piv = m.pivot_table(index="id", columns="partner", values="done",
                        aggfunc="first").fillna(False)
    return set(piv[(piv["He"]) & (piv["H2"])].index)


def build_pool(lib, man):
    ok = complete_ids(man)
    p = lib[lib["id"].isin(ok) & lib["B_cm1"].between(0.1, 1.0)].copy()
    Rg = np.linspace(3.5, 25.0, 6000)            # Bohr, outside inner barrier
    p["eps0_iso"] = p.apply(lambda r: -pg.v0(Rg, r.A, r.a, r.C6).min() * H, axis=1)
    p["r1"] = p["d_eps_odd"].abs() / p["eps0_iso"]
    p["r2"] = p["d_eps_even"].abs() / p["eps0_iso"]
    return p


def select(pool):
    chosen, used = [], set()
    for case, mode, dt in CASES:
        c = pool[pool["eps0_iso"].between(dt * 0.6, dt * 1.5) &
                 ~pool["id"].isin(used)].copy()
        if mode == "weakweak":
            c["score"] = -(c["r1"] + c["r2"])            # most isotropic
        elif mode == "v1large":
            c["score"] = c["r1"] - c["r2"]               # V1 up, V2 down
        else:                                            # v2large
            c["score"] = c["r2"] - c["r1"]
        # tie toward target depth
        c["score"] -= 0.15 * (c["eps0_iso"] - dt).abs() / dt
        best = c.sort_values("score", ascending=False).iloc[0]
        used.add(int(best["id"]))
        chosen.append((case, mode, dt, best))
    return chosen


def deexc_rates(out, mu, B):
    """{(i,f): k_markers[25]} for de-excitation transitions above KMIN."""
    df = parse_file(out)
    deexc = df[df["j_F"] < df["j_I"]]
    res = {}
    for (i, f), sub in deexc.groupby(["j_I", "j_F"]):
        sub = sub.sort_values("E")
        Ekin = sub["E"].values - B * i * (i + 1)
        m = Ekin > 0
        if m.sum() < 3:
            continue
        k = np.array([ek.k_of_T(Ekin[m], sub["SIG"].values[m], T, mu)
                      for T in T_MARK])
        if np.isfinite(k).any() and np.nanmax(k) >= KMIN:
            res[(int(i), int(f))] = k
    return res


def fit_line(k):
    good = k > 0
    coef = np.polyfit(T_MARK[good] ** (-1.0 / 6.0), np.log10(k[good]), FIT_DEG)
    return 10.0 ** np.polyval(coef, T_FINE ** (-1.0 / 6.0))


def main():
    lib = pd.read_csv(os.path.join(REPO, "potentials_library.csv"))
    man = pd.read_csv(os.path.join(PROD, "manifest.csv"))
    pool = build_pool(lib, man)
    chosen = select(pool)

    # selection table
    rows = []
    for case, mode, dt, b in chosen:
        rows.append({"case": case, "id": int(b["id"]),
                     "V0_depth_cm1": round(b["eps0_iso"], 1),
                     "B_cm1": round(b["B_cm1"], 3),
                     "V1/V0": round(b["r1"], 2), "V2/V0": round(b["r2"], 2),
                     "target_depth": dt, "anisotropy": mode})
    seltab = pd.DataFrame(rows)
    seltab.to_csv(os.path.join(HERE, "selected_cases.csv"), index=False)
    print(seltab.to_string(index=False))

    desc = {"weakweak": "V1 & V2 weak", "v1large": "V1 large, V2 weak",
            "v2large": "V1 weak, V2 large"}
    with PdfPages(PDF) as pdf:
        # summary page
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.axis("off")
        ax.set_title("Selected sample potentials (B = 0.1-1 cm$^{-1}$, completed runs)",
                     fontsize=13, weight="bold")
        tb = ax.table(cellText=seltab.values, colLabels=seltab.columns,
                      loc="center", cellLoc="center")
        tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 1.6)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        for case, mode, dt, b in chosen:
            pid = int(b["id"])
            decks = {}
            for pn, mucol in (("He", "mu_He"), ("H2", "mu_H2")):
                out = os.path.join(PROD, man[(man["id"] == pid) &
                                   (man["partner"] == pn)].iloc[0]["dir"],
                                   f"pot{pid:04d}_{pn}.out")
                decks[pn] = (deexc_rates(out, float(b[mucol]), float(b["B_cm1"])))
            # fixed low-j transitions (present in either partner); fall back to
            # strongest by peak if too few of the fixed set survived KMIN
            present = set(decks["He"]) | set(decks["H2"])
            top = [tr for tr in FIXED_TRANS if tr in present]
            if len(top) < 3:
                peak = {}
                for pn in decks:
                    for tr, k in decks[pn].items():
                        peak[tr] = max(peak.get(tr, 0), np.nanmax(k))
                extra = [t for t, _ in sorted(peak.items(), key=lambda kv: -kv[1])
                         if t not in top]
                top += extra[:5 - len(top)]

            fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharex=True)
            colors = plt.cm.tab10(np.linspace(0, 1, len(top)))
            for ax, pn in zip(axes, ("He", "H2")):
                for tr, col in zip(top, colors):
                    k = decks[pn].get(tr)
                    if k is None:
                        continue
                    ax.plot(T_FINE, fit_line(k), "-", color=col, lw=1.6)
                    ax.plot(T_MARK, k, "o", color=col, ms=4,
                            label=f"{tr[0]}$\\to${tr[1]}")
                ax.set_xlabel("T (K)")
                ax.set_ylabel("k  (cm$^3$/s)")
                ax.set_title(f"{pn}")
                ax.set_xlim(0, 50)
                ax.set_ylim(bottom=0)
                ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
                ax.legend(fontsize=8, title="i$\\to$f")
                ax.grid(alpha=0.3)
            fig.suptitle(f"Case {case}: {desc[mode]}  (V0~{dt} cm$^{{-1}}$)   "
                         f"pot{pid:04d}:  V0={b['eps0_iso']:.1f} cm$^{{-1}}$,  "
                         f"B={b['B_cm1']:.3f} cm$^{{-1}}$,  "
                         f"V1/V0={b['r1']:.2f},  V2/V0={b['r2']:.2f}",
                         fontsize=11, weight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            pdf.savefig(fig); plt.close(fig)
            print(f"  page done: case {case} pot{pid:04d} "
                  f"({len(top)} transitions)")

    print(f"\nPDF -> {PDF}")


if __name__ == "__main__":
    main()
