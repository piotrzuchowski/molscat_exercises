#!/usr/bin/env python3
"""Opacity functions sigma^J(i->f) vs J from MOLSCAT partial cross sections.

MOLSCAT prints (when ISIGPR>0 and IPRINT>=3) a block

  ** ... STATE-TO-STATE PARTIAL CROSS SECTIONS (ANGSTROM**2) FROM LEVEL I TO LEVEL F ...
        FOR JTOT = <J> AND SYMMETRY BLOCK = <IB> AT ENERGY(<n>) =<E> CM-1
     F  I =     <i1>        <i2>   ...        (column = source level I)
     <f1>     <sig(i1->f1)> <sig(i2->f1)> ... (row    = target level F)
     <f2>     ...

The cell at (row F, col I) is sigma FROM level I TO level F for that single
(JTOT, symmetry block).  The opacity function is the partial-wave contribution
per JTOT, so we SUM over all symmetry blocks sharing the same (JTOT, energy).
Level index L corresponds to rotor j = L-1 (JLEVEL = 0,1,2,...).

Output: opacity_functions.pdf  (one row of panels per potential; columns =
transitions 1->0, 2->1, 2->0; one curve per energy).
"""
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

# (lib id, short description) for the three runs
POTS = [("opac457_He.out", 457, "weak/weak (eps0~9)"),
        ("opac1466_He.out", 1466, "V1-large (eps0~106)"),
        ("opac1244_He.out", 1244, "V2-large (eps0~104)")]

# transitions as (j_i, j_f); level index = j+1
TRANS = [(1, 0), (2, 1), (2, 0)]

# physical constants (atomic units) for the centrifugal-barrier cutoff
HARTREE2CM = 219474.6313702
AMU2AU = 1822.888486209          # 1 amu in electron masses

# per-potential (C6 [a.u.], B [cm^-1], reduced mass mu_He [amu]) from the library
_lib = pd.read_csv(os.path.join(REPO, "potentials_library.csv")).set_index("id")
LIBINFO = {pid: (float(_lib.loc[pid, "C6"]), float(_lib.loc[pid, "B_cm1"]),
                 float(_lib.loc[pid, "mu_He"]))
           for _, pid, _ in POTS}


def j_barrier(E_cm1, ji, C6_au, B_cm1, mu_amu):
    """Partial wave whose van der Waals centrifugal barrier height equals the
    entrance-channel (level ji) collision energy.  Returns None if closed.

        Vmax^L = [hbar^2 L(L+1)/mu]^{3/2} (54 C6)^{-1/2}  = E_coll
      => L(L+1) = (mu/hbar^2) E_coll^{2/3} (54 C6)^{1/3}   (atomic units, hbar=1)
    """
    E_coll = E_cm1 - B_cm1 * ji * (ji + 1)      # entrance-channel kinetic energy
    if E_coll <= 0:
        return None
    E_au = E_coll / HARTREE2CM
    mu_au = mu_amu * AMU2AU
    LL1 = mu_au * E_au ** (2.0 / 3.0) * (54.0 * C6_au) ** (1.0 / 3.0)
    return 0.5 * (-1.0 + np.sqrt(1.0 + 4.0 * LL1))

HDR = re.compile(
    r"FOR JTOT\s*=\s*(\d+)\s+AND SYMMETRY BLOCK\s*=\s*(-?\d+).*?"
    r"AT ENERGY\(\s*(\d+)\)\s*=\s*([0-9.eE+-]+)")
BLOCK_MARK = "STATE-TO-STATE PARTIAL CROSS SECTIONS"


def parse(path):
    """-> dict[(Ecm1, jtot)][(j_i, j_f)] = sigma (summed over blocks)."""
    with open(path) as fh:
        lines = fh.readlines()
    # data[(inrg, E)][jtot][(i_lev, f_lev)] = sigma  (block-summed)
    acc = {}
    n = len(lines)
    k = 0
    while k < n:
        if BLOCK_MARK in lines[k]:
            # header with JTOT / block / energy is on the next 1-2 lines
            m = None
            for j in range(k, min(k + 4, n)):
                m = HDR.search(lines[j])
                if m:
                    break
            if not m:
                k += 1
                continue
            jtot = int(m.group(1))
            inrg = int(m.group(3))
            E = float(m.group(4))
            # find the column-header line "  F  I = ..."
            h = j + 1
            while h < n and "I =" not in lines[h]:
                h += 1
            if h >= n:
                k = j + 1
                continue
            cols = [int(x) for x in re.findall(r"\d+", lines[h].split("I =")[1])]
            # read matrix rows until a blank / separator line
            r = h + 1
            cell = {}
            while r < n:
                s = lines[r].rstrip("\n")
                if not s.strip() or s.strip().startswith("*"):
                    break
                toks = s.split()
                # row label (F level) + values, all numeric
                if not re.match(r"^\s*\d+\s", s):
                    break
                try:
                    flev = int(toks[0])
                    vals = [float(t) for t in toks[1:1 + len(cols)]]
                except ValueError:
                    break
                if len(vals) != len(cols):
                    break
                for ci, ilev in enumerate(cols):
                    cell[(ilev, flev)] = vals[ci]   # sigma from I=ilev to F=flev
                r += 1
            key = (inrg, E)
            d = acc.setdefault(key, {}).setdefault(jtot, {})
            for (ilev, flev), v in cell.items():
                d[(ilev, flev)] = d.get((ilev, flev), 0.0) + v
            k = r
        else:
            k += 1
    # collapse to dict[(E, jtot)][(j_i, j_f)]
    out = {}
    for (inrg, E), perj in acc.items():
        for jtot, cell in perj.items():
            tgt = out.setdefault((E, jtot), {})
            for (ilev, flev), v in cell.items():
                tgt[(ilev - 1, flev - 1)] = v       # level -> j
    return out


def opacity(parsed, ji, jf):
    """-> dict[E] = (Jarr, sigarr) for transition ji->jf."""
    byE = {}
    for (E, jtot), cell in parsed.items():
        if (ji, jf) in cell:
            byE.setdefault(E, []).append((jtot, cell[(ji, jf)]))
    out = {}
    for E, pts in byE.items():
        pts.sort()
        J = np.array([p[0] for p in pts])
        sg = np.array([p[1] for p in pts])
        out[E] = (J, sg)
    return out


def main():
    parsed_all = {}
    for fn, pid, desc in POTS:
        p = os.path.join(HERE, fn)
        if not os.path.exists(p):
            print(f"MISSING {fn}")
            continue
        parsed_all[pid] = (parse(p), desc)
        es = sorted({E for (E, _) in parsed_all[pid][0]})
        print(f"pot{pid}: {len(es)} energies, "
              f"{len(parsed_all[pid][0])} (E,JTOT) blocks; "
              f"E={[round(e,2) for e in es]}")

    pdf_path = os.path.join(HERE, "opacity_functions.pdf")
    with PdfPages(pdf_path) as pdf:
        for pid, (parsed, desc) in parsed_all.items():
            energies = sorted({E for (E, _) in parsed})
            cmap = plt.cm.viridis(np.linspace(0, 1, len(energies)))
            fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
            C6_au, B_cm1, mu_amu = LIBINFO[pid]
            for ax, (ji, jf) in zip(axes, TRANS):
                op = opacity(parsed, ji, jf)
                peak = 0.0
                for E, col in zip(energies, cmap):
                    if E not in op:
                        continue
                    J, sg = op[E]
                    pos = sg > 0                      # log scale: drop non-positive
                    peak = max(peak, sg[pos].max() if pos.any() else 0.0)
                    ax.plot(J[pos], sg[pos], "-", color=col, lw=1.2, marker="o",
                            ms=2.5, label=f"{E:.2f}")
                    # vertical line: J where vdW centrifugal barrier == E_coll
                    jc = j_barrier(E, ji, C6_au, B_cm1, mu_amu)
                    if jc is not None:
                        ax.axvline(jc, color=col, ls="--", lw=1.0, alpha=0.7)
                ax.set_xlabel("J (total angular momentum)")
                ax.set_ylabel(r"$\sigma^J$ ($\mathrm{\AA}^2$)")
                ax.set_title(f"{ji}$\\to${jf}")
                ax.set_yscale("log")
                ax.set_ylim(bottom=1e-6, top=(peak * 3.0 if peak > 0 else None))
                ax.set_xlim(0, 40)
                ax.grid(alpha=0.3, which="both")
            axes[-1].legend(fontsize=6, title="E (cm$^{-1}$)", ncol=2,
                            loc="upper right")
            fig.suptitle(f"Opacity functions  pot{pid}  [{desc}]  (He partner)",
                         fontsize=12, weight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            pdf.savefig(fig)
            plt.close(fig)
            print(f"  page done: pot{pid}")
    print(f"\nPDF -> {pdf_path}")


if __name__ == "__main__":
    main()
