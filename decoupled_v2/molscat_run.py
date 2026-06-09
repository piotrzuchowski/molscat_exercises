#!/usr/bin/env python3
"""Drive MOLSCAT (ITYPE=1, atom + linear rigid rotor) over the potential library.

For each library row we build a MOLSCAT 2020 namelist deck in which the Born-Mayer
+ dispersion potential is supplied analytically:
    lambda=0 : A*exp(-a R) - C6 R^-6
    lambda=1 : eta*A*exp(-a R)  - a1*C6 R^-6   (odd:  wall eta, tail a1)
    lambda=2 : zeta*A*exp(-a R) - a2*C6 R^-6   (even: wall zeta, tail a2)
The exponential term uses NPOWER=0 with exponent E=-a; the dispersion term uses
NPOWER=-6. Working units are cm^-1 / Angstrom, so the au parameters are converted.

Per case we compute state-to-state integral cross sections at a few collision
energies scaled to the rotor spacing B (so the open-channel structure is the same
for every system), with the H2 reduced mass. Validation default: first 5 rows.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

HARTREE2CM = 219474.6313702
BOHR2ANG = 0.529177210903
HBAR2_2MU = 16.85762920            # cm^-1 * amu * Ang^2 : E = HBAR2_2MU * k^2 / mu

CONFIG = {
    "binary": "/Users/piotrzuchowski/develop/molscat/source_code/molscat-basic",
    "library": "potentials_library.csv",
    "workdir": "molscat_runs",
    "results_csv": "molscat_xsec.csv",
    "n_cases": 5,                  # validation batch; raise to 200 once happy
    "mu_col": "mu_H2",            # H2 partner
    "energy_factors": (4.0, 16.0, 64.0),   # collision energies = factor * B  [cm^-1]
    "RMIN_bohr": 4.0,             # generator-screened outer edge; RMIN never exceeds it
    "RMAX_bohr": 30.0,
    # Adaptive RMIN guard: start propagation deep in the classically forbidden
    # region so the BCYOMN hard wall is exact and can never mask a genuinely open
    # channel. Require V0(RMIN) >= wall_safety * E_max (worst case is J=0,L=0).
    "wall_safety": 30.0,          # V0(RMIN) must exceed this multiple of E_max
    "RMIN_floor_bohr": 2.5,       # hard inner limit; never start inside this
    "n_scan": 400,                # grid for the no-hole re-check on [rmin, RMIN_bohr]
    "DR": 0.01,                   # propagator step [Ang]
    "IPROPS": 6,
    "JMAX": 10,                   # rotor basis (open up to J~7 at 64B, + buffer)
    "ltop_buffer": 12,            # extra partial waves above the semiclassical l_max
}


# --------------------------------------------------------------------------- #
def au_to_molscat(A_au, a_au, C6_au):
    """Convert (A [Hartree], a [Bohr^-1], C6 [au]) to cm^-1 / Angstrom units."""
    A_cm = A_au * HARTREE2CM
    a_ang = a_au / BOHR2ANG                  # exp(-a R): coeff in Ang^-1
    C6_cm = C6_au * HARTREE2CM * BOHR2ANG ** 6
    return A_cm, a_ang, C6_cm


def estimate_jtotu(E_max, mu, rmax_ang):
    """Semiclassical l_max = R_max * k, k = sqrt(E*mu/HBAR2_2MU), plus a buffer."""
    k = np.sqrt(E_max * mu / HBAR2_2MU)
    return int(np.ceil(rmax_ang * k)) + CONFIG["ltop_buffer"]


def v0_cm(R_bohr, A_au, a_au, C6_au):
    """Isotropic Born-Mayer + dispersion potential [cm^-1], R in Bohr (au params)."""
    return (A_au * np.exp(-a_au * R_bohr) - C6_au * R_bohr ** -6) * HARTREE2CM


def dv0_dR_au(R_bohr, A_au, a_au, C6_au):
    """Radial derivative of v0; > 0 marks the inner turnover toward the collapse."""
    return -a_au * A_au * np.exp(-a_au * R_bohr) + 6.0 * C6_au * R_bohr ** -7


def choose_rmin(row, cfg):
    """Pick RMIN [Bohr] so V0(RMIN) >= wall_safety * E_max (deep forbidden region).

    Returns (rmin_bohr, None) on success, or (None, reason) if the repulsive wall
    is too short to support E_max before the unphysical -C6/R^6 collapse -- a hard
    failure rather than a silently-masked open channel.
    """
    A, a, C6 = row["A"], row["a"], row["C6"]
    e_max = max(cfg["energy_factors"]) * row["B_cm1"]
    target = cfg["wall_safety"] * e_max          # required V0(RMIN) [cm^-1]
    r_edge, r_floor = cfg["RMIN_bohr"], cfg["RMIN_floor_bohr"]

    if v0_cm(r_edge, A, a, C6) >= target:
        return r_edge, None                      # screened edge already deep enough
    if v0_cm(r_floor, A, a, C6) < target:
        return None, (f"wall_too_short: V0({r_floor:.2f}a0)="
                      f"{v0_cm(r_floor, A, a, C6):.1f} < {cfg['wall_safety']:.0f}"
                      f"*Emax={target:.1f} cm-1")
    rmin = brentq(lambda R: v0_cm(R, A, a, C6) - target, r_floor, r_edge)
    # Re-screen the inward extension [rmin, r_edge] the generator did not check:
    # V0 must stay positive and monotonically repulsive (no hole) down to rmin.
    scan = np.linspace(rmin, r_edge, cfg["n_scan"])
    if np.any(v0_cm(scan, A, a, C6) <= 0.0) or np.any(dv0_dR_au(scan, A, a, C6) > 0.0):
        return None, f"inner_collapse_below_{r_edge:.1f}a0 (rmin={rmin:.3f})"
    return rmin, None


def build_deck(row, cfg, rmin_bohr):
    """Return the MOLSCAT input deck (string) for one library row."""
    A_cm, a_ang, C6_cm = au_to_molscat(row["A"], row["a"], row["C6"])
    B = row["B_cm1"]
    mu = row[cfg["mu_col"]]
    energies = [f * B for f in cfg["energy_factors"]]
    jtotu = estimate_jtotu(max(energies), mu, cfg["RMAX_bohr"] * BOHR2ANG)

    a1, a2, eta = row["a1"], row["a2"], row["eta"]
    zeta = row["zeta"] if "zeta" in row else a2   # old libs: V2 = a2*V0 (zeta=a2)
    # Flattened term arrays: lambda = 0,1,2 ; each has [exp term, R^-6 term].
    # Wall (exp) carries (eta, zeta); dispersion tail (R^-6) carries (a1, a2).
    A_arr = [A_cm, -C6_cm, eta * A_cm, -a1 * C6_cm, zeta * A_cm, -a2 * C6_cm]
    E_arr = [-a_ang, -a_ang, -a_ang]            # one per NPOWER==0 entry

    en_str = ", ".join(f"{e:.6f}" for e in energies)
    A_str = ", ".join(f"{x:.8E}" for x in A_arr)
    E_str = ", ".join(f"{x:.8E}" for x in E_arr)

    # BCYOMN: hard-wall log-derivative for channels MOLSCAT flags "locally open"
    # at RMIN. With odd-J rotor states the lambda=0 coupling carries a (-1)^J sign,
    # so odd-J diagonals go negative under a repulsive V0 and trip the WKB-closed
    # initialiser (IRMSET=0). At RMIN the wavefunction is ~0 (V0 >> E), so a
    # hard wall is exact -- verified: sigma(0->0) matches the single-channel run.
    return f""" &INPUT
  LABEL  = 'pot id={int(row["id"])} ({cfg["mu_col"]})',
  URED   = {mu:.6f},
  IPRINT = 1, ISIGPR = 2,
  RMIN   = {rmin_bohr * BOHR2ANG:.5f}, RMAX = {cfg["RMAX_bohr"] * BOHR2ANG:.5f},
  IRMSET = 0, BCYOMN = 1.0E6,
  IPROPS = {cfg["IPROPS"]}, DR = {cfg["DR"]},
  NNRG   = {len(energies)}, ENERGY = {en_str},
  JTOTL  = 0, JTOTU = {jtotu}, JSTEP = 1,
 /

 &BASIS
  ITYPE  = 1,
  ROTI   = {B:.8f},
  JMAX   = {cfg["JMAX"]},
 /

 &POTL
  MXLAM  = 3, LAMBDA = 0, 1, 2,
  NTERM  = 2, 2, 2,
  NPOWER = 0, -6, 0, -6, 0, -6,
  E      = {E_str},
  A      = {A_str},
  RM     = 1.0, EPSIL = 1.0,
 /
"""


# --------------------------------------------------------------------------- #
_SIG_HDR = re.compile(r"SIG\(F,I\)")


def parse_sigma(text):
    """Parse state-to-state integral cross sections from MOLSCAT stdout.

    Returns a list of dicts: energy, F, I, J_F, J_I, sigma [Ang^2].
    """
    rows, in_block = [], False
    for line in text.splitlines():
        if _SIG_HDR.search(line):
            in_block = True
            continue
        if in_block:
            if "TOTAL INELASTIC" in line or "MOLSCAT" in line:
                in_block = False
                continue
            tok = line.split()
            if len(tok) == 7:
                try:
                    energy = float(tok[0])
                    f, i = int(tok[4]), int(tok[5])
                    sig = float(tok[6])
                except ValueError:
                    continue
                rows.append({"energy": energy, "F": f, "I": i,
                             "J_F": f - 1, "J_I": i - 1, "sigma": sig})
    return rows


# --------------------------------------------------------------------------- #
def run_case(row, cfg, wd):
    rmin_bohr, err = choose_rmin(row, cfg)
    if err:                                  # wall can't support E_max: do not run
        return None, None, err
    deck = build_deck(row, cfg, rmin_bohr)
    cid = int(row["id"])
    inp = wd / f"pot_{cid:05d}.inp"
    out = wd / f"pot_{cid:05d}.out"
    inp.write_text(deck)
    proc = subprocess.run([cfg["binary"]], stdin=open(inp),
                          capture_output=True, text=True)
    out.write_text(proc.stdout + "\n--- STDERR ---\n" + proc.stderr)
    if proc.returncode != 0:
        return None, rmin_bohr, f"returncode={proc.returncode}"
    sig = parse_sigma(proc.stdout)
    if not sig:
        return None, rmin_bohr, "no_cross_sections_parsed"
    return sig, rmin_bohr, None


def main(cfg=CONFIG):
    lib = pd.read_csv(cfg["library"]).head(cfg["n_cases"])
    wd = Path(cfg["workdir"])
    wd.mkdir(exist_ok=True)

    all_rows, status = [], []
    for _, row in lib.iterrows():
        cid = int(row["id"])
        sig, rmin_bohr, err = run_case(row, cfg, wd)
        if err:
            status.append((cid, "FAIL", err, 0))
            print(f"  id={cid:5d}  FAIL  {err}")
            continue
        for s in sig:
            s["id"] = cid
            s["B_cm1"] = row["B_cm1"]
            s["mu"] = row[cfg["mu_col"]]
            s["kind"] = "elastic" if s["F"] == s["I"] else "inelastic"
            all_rows.append(s)
        n_inel = sum(1 for s in sig if s["F"] != s["I"])
        status.append((cid, "ok", "", len(sig)))
        flag = "" if rmin_bohr == cfg["RMIN_bohr"] else f", RMIN->{rmin_bohr:.3f}a0"
        print(f"  id={cid:5d}  ok    {len(sig):3d} cross sections "
              f"({n_inel} inelastic), B={row['B_cm1']:.4g} cm-1{flag}")

    if all_rows:
        df = pd.DataFrame(all_rows)[
            ["id", "B_cm1", "mu", "energy", "I", "F", "J_I", "J_F", "kind", "sigma"]]
        df.to_csv(cfg["results_csv"], index=False)
        n_ok = sum(1 for s in status if s[1] == "ok")
        print(f"\n{n_ok}/{len(status)} cases succeeded; "
              f"{len(df)} cross sections -> {cfg['results_csv']}")
        print(f"Decks + raw output in {cfg['workdir']}/")
    else:
        print("\nNo cross sections produced; inspect decks in", cfg["workdir"])
    return all_rows, status


if __name__ == "__main__":
    if len(sys.argv) > 1:
        CONFIG["n_cases"] = int(sys.argv[1])
    main()
