#!/usr/bin/env python3
"""Build MOLSCAT ITYPE=1 input decks from library potentials (smoke test).

For each selected potential and each collision partner (He, H2) writes a
molscat-basic-style input deck using:
  * the analytic &POTL form (exact map of Born-Mayer + dispersion, no recompile)
  * dynamic RMIN  (radius where the repulsive wall reaches E_REF, all angles)
  * dynamic RMAX  (radius where the isotropic tail drops below TOL_CM)
  * Jmax from jmax_calc.set_jmax(B, Emax)
  * energy grid: 30 log points 0.1-50 cm^-1, then even continuation with the
    last log-step up to E_MAX, passed as an explicit ENERGY() list (DNRG=0)

Potential -> &POTL mapping (RM=1 Angstrom, EPSIL=1 cm^-1, so R in Ang, V in cm^-1):
  lambda=0: A*exp(-a R) - C6 R^-6
  lambda=1: eta*A*exp(-a R) - a1*C6 R^-6
  lambda=2: a2*A*exp(-a R) - a2*C6 R^-6
with NPOWER=0 -> exponential term A*exp(E*R), NPOWER=-6 -> power term A*R^-6.
Internal generator units are atomic (R in Bohr, energy Hartree), converted here.
"""

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import potential_generator as pg
import jmax_calc

HARTREE2CM = pg.HARTREE2CM
BOHR2ANG = pg.BOHR2ANG

# ----------------------------- configuration ------------------------------ #
LIBRARY = "../potentials_library.csv"
SELECT_IDS = [39, 79, 82, 94, 96]      # 5 tractable potentials (B>=4 -> jmax<=12)
PARTNERS = {"He": "mu_He", "H2": "mu_H2"}

E_LO, E_HI_LOG, N_LOG = 0.1, 50.0, 30   # log grid
E_MAX = 200.0                            # even-continuation upper bound
E_REF_RMIN = 1.0e4                       # cm^-1, wall energy defining RMIN
TOL_RMAX = 1.0e-3                        # cm^-1, tail cutoff defining RMAX
RMAX_LO, RMAX_HI = 30.0, 80.0            # Angstrom clamp
RMAX_FIXED = None                        # if set, override dynamic RMAX [Ang]
DR = 0.01                                # propagator step (Angstrom)
JMAX_CAP = 30
JSTEP = 1                                # partial-wave step (override via argv)
DTOL = None                              # elastic conv tol [Ang^2]; None->molscat default 0.3
OUTDIR = "."                             # where decks are written (override via argv)


# ------------------------------ energy grid -------------------------------- #
def energy_grid():
    log = np.logspace(np.log10(E_LO), np.log10(E_HI_LOG), N_LOG)
    step = log[-1] - log[-2]
    cont, e = [], log[-1] + step
    while e <= E_MAX + 1e-9:
        cont.append(e)
        e += step
    return np.concatenate([log, cont])


# --------------------------- dynamic RMIN / RMAX --------------------------- #
def dynamic_rmin(p):
    """Min over angles of the radius where V(R,theta) drops through E_REF_RMIN."""
    eref_h = E_REF_RMIN / HARTREE2CM
    Rg = np.linspace(1.0, 15.0, 8000)           # Bohr
    crossings = []
    for ct in (1.0, 0.0, -1.0):
        V = pg.v_theta(Rg, p["A"], p["a"], p["C6"], p["a1"], p["a2"], p["eta"], ct)
        idx = np.where(V >= eref_h)[0]
        if idx.size:
            crossings.append(Rg[idx[-1]])
    rmin_bohr = min(crossings) if crossings else 1.5
    return rmin_bohr * BOHR2ANG


def dynamic_rmax(p):
    """Radius (Ang) where the isotropic -C6/R^6 tail falls below TOL_RMAX."""
    if RMAX_FIXED is not None:
        return float(RMAX_FIXED)
    c6term = HARTREE2CM * p["C6"] * BOHR2ANG ** 6    # cm^-1 * Ang^6
    R = (c6term / TOL_RMAX) ** (1.0 / 6.0)
    return float(np.clip(R, RMAX_LO, RMAX_HI))


# ------------------------------ &POTL mapping ------------------------------ #
def potl_coeffs(p):
    """Return (Eexp, [A0e,A0p,A1e,A1p,A2e,A2p]) for the analytic &POTL block."""
    Eexp = -p["a"] / BOHR2ANG                         # Ang^-1 (exp argument)
    c6term = HARTREE2CM * p["C6"] * BOHR2ANG ** 6     # cm^-1 Ang^6
    Apre = HARTREE2CM * p["A"]                        # cm^-1 (exp prefactor)
    A0e, A0p = Apre, -c6term
    A1e, A1p = p["eta"] * Apre, -p["a1"] * c6term
    A2e, A2p = p["a2"] * Apre, -p["a2"] * c6term
    return Eexp, [A0e, A0p, A1e, A1p, A2e, A2p]


# ------------------------------ deck writer -------------------------------- #
def fmt_list(vals, fmt, per_line=6, indent=13):
    out, pad = [], " " * indent
    for i in range(0, len(vals), per_line):
        chunk = ", ".join(fmt % v for v in vals[i:i + per_line])
        out.append(chunk)
    return (",\n" + pad).join(out)


def build_deck(p, partner, ured, egrid, jmax):
    rmin = dynamic_rmin(p)
    rmax = dynamic_rmax(p)
    Eexp, A = potl_coeffs(p)
    jlist = list(range(jmax + 1))
    nlevel = len(jlist)

    elist = fmt_list(list(egrid), "%.6f", per_line=6, indent=29)
    jstr = ", ".join("%d" % j for j in jlist)
    astr = fmt_list(A, "%.8E", per_line=2, indent=20)
    dtol_str = f"   DTOL = {DTOL}," if DTOL is not None else ""

    deck = f""" &INPUT
    LABEL  = 'lib id {int(p['id'])}  partner {partner}  (smoke test)',
    URED   = {ured:.6f},
    IPRINT = 1,    ISIGPR = 2,
    RMIN   = {rmin:.4f},  RMAX   = {rmax:.4f},  IRMSET = 0,
    IPROPS = 6,    DR     = {DR},
    JTOTL  = 0,    JTOTU  = 99999,   JSTEP = {JSTEP},{dtol_str}
    NNRG   = {len(egrid)},
    ENERGY = {elist},
 /

 &BASIS
    ITYPE  = 1,   BE     = {p['B_cm1']:.8f},
    NLEVEL = {nlevel},
    JLEVEL = {jstr},
 /

 &POTL
    MXLAM  = 3,   LAMBDA = 0, 1, 2,
    NTERM  = 2, 2, 2,
    NPOWER = 0, -6,  0, -6,  0, -6,
    A      = {astr},
    E      = {Eexp:.8E}, {Eexp:.8E}, {Eexp:.8E},
    RM     = 1.0, EPSIL  = 1.0,
 /
"""
    return deck, dict(id=int(p["id"]), partner=partner, ured=ured,
                      rmin=rmin, rmax=rmax, jmax=jmax, nlevel=nlevel,
                      B=p["B_cm1"], nnrg=len(egrid))


def main():
    global JSTEP, DTOL, OUTDIR, DR, RMAX_FIXED
    here = os.path.dirname(os.path.abspath(__file__))
    # CLI: make_decks.py [outdir] [jstep] [dtol] [dr] [rmax]
    if len(sys.argv) > 1:
        OUTDIR = sys.argv[1]
    if len(sys.argv) > 2:
        JSTEP = int(sys.argv[2])
    if len(sys.argv) > 3:
        DTOL = float(sys.argv[3])
    if len(sys.argv) > 4:
        DR = float(sys.argv[4])
    if len(sys.argv) > 5:
        RMAX_FIXED = float(sys.argv[5])
    libpath = os.path.join(here, LIBRARY)
    outdir = os.path.join(here, OUTDIR)
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)
    print(f"Settings: OUTDIR={OUTDIR}  JSTEP={JSTEP}  DTOL={DTOL}  "
          f"DR={DR}  RMAX_FIXED={RMAX_FIXED}")
    df = pd.read_csv(libpath)
    egrid = energy_grid()
    jmax_E = float(egrid.max())
    print(f"Energy grid: {len(egrid)} points, "
          f"{egrid[0]:.3f} .. {egrid[-1]:.3f} cm^-1 "
          f"(log step continuation {egrid[-1]-egrid[-2]:.3f}); "
          f"jmax reference E = {jmax_E:.1f}")

    manifest = []
    for pid in SELECT_IDS:
        row = df.loc[df["id"] == pid]
        if row.empty:
            print(f"  !! id {pid} not in library, skipping")
            continue
        p = row.iloc[0]
        jmax = jmax_calc.set_jmax(p["B_cm1"], jmax_E, cap=JMAX_CAP)
        for partner, mucol in PARTNERS.items():
            deck, meta = build_deck(p, partner, float(p[mucol]), egrid, jmax)
            fname = f"pot{pid:04d}_{partner}.input"
            with open(fname, "w") as fh:
                fh.write(deck)
            manifest.append({**meta, "input": fname,
                             "output": fname.replace(".input", ".out")})
            print(f"  wrote {fname}: jmax={jmax} nlevel={meta['nlevel']} "
                  f"RMIN={meta['rmin']:.3f} RMAX={meta['rmax']:.2f} URED={meta['ured']:.4f}")

    pd.DataFrame(manifest).to_csv("manifest.csv", index=False)
    print(f"\nWrote {len(manifest)} decks + manifest.csv")


if __name__ == "__main__":
    main()
