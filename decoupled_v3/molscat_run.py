#!/usr/bin/env python3
"""decoupled_v3 MOLSCAT deck helpers (ITYPE=1, atom + linear rotor).

lambda decomposition supplied analytically (cm^-1 / Angstrom):
    lambda=0 :      A e^{-aR}  -      C6 R^-6        (NPOWER  0, -6)
    lambda=1 : s1*( A1 e^{-b1 R} - C7 R^-7 )         (NPOWER  0, -7)  <-- odd, R^-7
    lambda=2 : zeta A e^{-aR}  - a2   C6 R^-6        (NPOWER  0, -6)
The lambda=1 exponential uses its OWN exponent b1 (independent of a). s1=+-1 is the
internal odd-anisotropy sign carried through from the library.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq

HARTREE2CM = 219474.6313702
BOHR2ANG = 0.529177210903
HBAR2_2MU = 16.85762920

CONFIG = {
    "binary": "/home/pzuch007/molscat_exercises/source_code/molscat-basic",
    "library": "potentials_library.csv",
    "mu_col": "mu_H2",
    "RMIN_bohr": 4.0, "RMAX_bohr": 30.0,
    "wall_safety": 30.0, "RMIN_floor_bohr": 2.5, "n_scan": 400,
    "DR": 0.04, "DTOL": 10.0, "IPROPS": 6, "JMAX": 15, "ltop_buffer": 12,
}


def au_to_molscat_exp(A_au, exp_au):
    """A [Hartree]->cm^-1 prefactor; exponent [1/Bohr]->1/Ang."""
    return A_au * HARTREE2CM, exp_au / BOHR2ANG


def disp_to_molscat(Cn_au, n):
    """Dispersion coeff Cn [Hartree*Bohr^n] -> cm^-1 * Ang^n for R^-n term."""
    return Cn_au * HARTREE2CM * BOHR2ANG ** n


def estimate_jtotu(E_max, mu, rmax_ang):
    k = np.sqrt(E_max * mu / HBAR2_2MU)
    return int(np.ceil(rmax_ang * k)) + CONFIG["ltop_buffer"]


def v0_cm(R_bohr, A, a, C6):
    return (A * np.exp(-a * R_bohr) - C6 * R_bohr ** -6) * HARTREE2CM


def dv0_dR_au(R_bohr, A, a, C6):
    return -a * A * np.exp(-a * R_bohr) + 6.0 * C6 * R_bohr ** -7


def choose_rmin(row, cfg):
    """RMIN [Bohr] so isotropic V0(RMIN) >= wall_safety*E_max. BCYOMN handles any
    odd-anisotropy 'locally open' channels at RMIN (hard wall)."""
    A, a, C6 = row["A"], row["a"], row["C6"]
    e_max = cfg["_emax_cm1"]
    target = cfg["wall_safety"] * e_max
    r_edge, r_floor = cfg["RMIN_bohr"], cfg["RMIN_floor_bohr"]
    if v0_cm(r_edge, A, a, C6) >= target:
        return r_edge, None
    if v0_cm(r_floor, A, a, C6) < target:
        return None, (f"wall_too_short: V0({r_floor:.2f})="
                      f"{v0_cm(r_floor, A, a, C6):.1f} < {target:.1f}")
    rmin = brentq(lambda R: v0_cm(R, A, a, C6) - target, r_floor, r_edge)
    scan = np.linspace(rmin, r_edge, cfg["n_scan"])
    if np.any(v0_cm(scan, A, a, C6) <= 0.0) or np.any(dv0_dR_au(scan, A, a, C6) > 0.0):
        return None, f"inner_collapse_below_{r_edge:.1f}a0"
    return rmin, None


def build_deck(row, cfg, rmin_bohr, energies):
    """MOLSCAT deck (string) with explicit dense ENERGY list for k(T)."""
    B = row["B_cm1"]
    mu = row[cfg["mu_col"]]
    jtotu = estimate_jtotu(max(energies), mu, cfg["RMAX_bohr"] * BOHR2ANG)

    A_cm, a_ang = au_to_molscat_exp(row["A"], row["a"])
    C6_cm = disp_to_molscat(row["C6"], 6)
    a2, zeta, s1 = row["a2"], row["zeta"], row["s1"]
    # lambda=1 independent term: s1*(A1 e^{-b1 R} - C7 R^-7)
    A1_cm, b1_ang = au_to_molscat_exp(row["A1_au"], row["b1_au"])
    C7_cm = disp_to_molscat(row["C7_au"], 7)

    # flattened term arrays, lambda = 0,1,2
    A_arr = [A_cm, -C6_cm,
             s1 * A1_cm, -s1 * C7_cm,
             zeta * A_cm, -a2 * C6_cm]
    E_arr = [-a_ang, -b1_ang, -a_ang]          # one exponent per NPOWER==0 term
    en_str = ", ".join(f"{e:.6f}" for e in energies)
    A_str = ", ".join(f"{x:.8E}" for x in A_arr)
    E_str = ", ".join(f"{x:.8E}" for x in E_arr)

    return f""" &INPUT
  LABEL  = 'pot id={int(row["id"])} ({cfg["mu_col"]}) v3 k(T)',
  URED   = {mu:.6f},
  IPRINT = 1, ISIGPR = 2,
  RMIN   = {rmin_bohr * BOHR2ANG:.5f}, RMAX = {cfg["RMAX_bohr"] * BOHR2ANG:.5f},
  IRMSET = 0, BCYOMN = 1.0E6,
  IPROPS = {cfg["IPROPS"]}, DR = {cfg["DR"]},
  NNRG   = {len(energies)}, ENERGY = {en_str},
  JTOTL  = 0, JTOTU = {jtotu}, JSTEP = 1, DTOL = {cfg["DTOL"]},
 /

 &BASIS
  ITYPE  = 1,
  ROTI   = {B:.8f},
  JMAX   = {cfg["JMAX"]},
 /

 &POTL
  MXLAM  = 3, LAMBDA = 0, 1, 2,
  NTERM  = 2, 2, 2,
  NPOWER = 0, -6, 0, -7, 0, -6,
  E      = {E_str},
  A      = {A_str},
  RM     = 1.0, EPSIL = 1.0,
 /
"""


# state-to-state integral cross-section parser (same format as v2)
import re
_SIG_HDR = re.compile(r"SIG\(F,I\)")


def parse_sigma(text):
    rows, inblk = [], False
    for line in text.splitlines():
        if _SIG_HDR.search(line):
            inblk = True
            continue
        if inblk:
            if "TOTAL INELASTIC" in line or "MOLSCAT" in line:
                inblk = False
                continue
            t = line.split()
            if len(t) == 7:
                try:
                    e = float(t[0]); f = int(t[4]); i = int(t[5]); s = float(t[6])
                except ValueError:
                    continue
                rows.append({"energy": e, "F": f, "I": i,
                             "J_F": f - 1, "J_I": i - 1, "sigma": s})
    return rows
