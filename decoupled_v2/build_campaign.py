#!/usr/bin/env python3
"""Build the decoupled_v2 cross-section campaign for the cluster.

Uses molscat_run.py's deck logic (decoupled V2 via zeta, BCYOMN hard wall,
adaptive RMIN, ROTI/JMAX) but writes one deck per accepted potential and a
SLURM array script + manifest, instead of running molscat in-process (the
cluster binary needs the MKL modules, so it must run inside a batch job).

Usage:
    python build_campaign.py [n_cases]      # default: all accepted potentials
Produces in production/:
    decks/pot_XXXXX.inp        one MOLSCAT deck per runnable potential
    inputs.list                ordered list of decks (array index -> deck)
    manifest.csv               id, rmin_bohr, B, mu, jtotu, n_energy, inp, out
    skipped.csv                potentials choose_rmin rejected (+ reason)
    run_campaign.sbatch        array job (throttled to <=100 concurrent)
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import molscat_run as mr           # deck + rmin logic (decoupled_v2)
import energy_kinteg_helper as ek  # dense k(T) energy grid (jmax15 machinery)

BINARY = "/home/pzuch007/molscat_exercises/source_code/molscat-basic"
ACCOUNT = "hpc-pzuch007-1692967268"
PARTITION = "lem-cpu-short"
MAX_CONCURRENT = 180
JMAX = 10                          # rotor basis (matches molscat_run CONFIG)
T_HI = 50.0                        # k(T) valid up to this T (E_max = 6 kB T_hi)
PROD = os.path.join(HERE, "production")
DECKS = os.path.join(PROD, "decks")


def avoid_thresholds(egrid, B, jmax, rel=2.0e-3):
    """Nudge energies sitting exactly on a rotor threshold B*j(j+1) off it, to
    dodge the RBESSK kappa=0 crash (an exactly-at-threshold channel)."""
    E = np.asarray(egrid, float).copy()
    thr = np.array([B * j * (j + 1) for j in range(1, jmax + 1)])
    for t in thr:
        if t <= 0:
            continue
        hit = np.abs(E - t) <= rel * t
        E[hit] = t * (1.0 - rel)
    return np.unique(E)


def dense_energies(B):
    """jmax15-style absolute collision-energy grid [cm^-1], threshold-safe."""
    eg = ek.energy_grid(B, JMAX, T_hi=T_HI)
    return avoid_thresholds(eg, B, JMAX)


def build_kt_deck(row, cfg, rmin_bohr, energies):
    """decoupled_v2 deck (zeta, BCYOMN, adaptive RMIN) with an explicit dense
    ENERGY list for k(T).  Mirrors molscat_run.build_deck otherwise."""
    A_cm, a_ang, C6_cm = mr.au_to_molscat(row["A"], row["a"], row["C6"])
    B = row["B_cm1"]
    mu = row[cfg["mu_col"]]
    jtotu = mr.estimate_jtotu(max(energies), mu, cfg["RMAX_bohr"] * mr.BOHR2ANG)
    a1, a2, eta = row["a1"], row["a2"], row["eta"]
    zeta = row["zeta"] if "zeta" in row else a2
    A_arr = [A_cm, -C6_cm, eta * A_cm, -a1 * C6_cm, zeta * A_cm, -a2 * C6_cm]
    en_str = ", ".join(f"{e:.6f}" for e in energies)
    A_str = ", ".join(f"{x:.8E}" for x in A_arr)
    E_str = ", ".join(f"{-a_ang:.8E}" for _ in range(3))
    return f""" &INPUT
  LABEL  = 'pot id={int(row["id"])} ({cfg["mu_col"]}) k(T)',
  URED   = {mu:.6f},
  IPRINT = 1, ISIGPR = 2,
  RMIN   = {rmin_bohr * mr.BOHR2ANG:.5f}, RMAX = {cfg["RMAX_bohr"] * mr.BOHR2ANG:.5f},
  IRMSET = 0, BCYOMN = 1.0E6,
  IPROPS = {cfg["IPROPS"]}, DR = {cfg["DR"]},
  NNRG   = {len(energies)}, ENERGY = {en_str},
  JTOTL  = 0, JTOTU = {jtotu}, JSTEP = 1, DTOL = {cfg.get("DTOL", 10.0)},
 /

 &BASIS
  ITYPE  = 1,
  ROTI   = {B:.8f},
  JMAX   = {JMAX},
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


def main():
    cfg = dict(mr.CONFIG)
    cfg["binary"] = BINARY
    # k(T) numerics: match the validated jmax15 production (tractable over the
    # dense 82-energy grid). DR=0.04 (4x fewer steps, lossless) + DTOL=10 (loose
    # elastic tol -> partial-wave sum terminates at convergence, not the cap).
    cfg["DR"] = 0.04
    cfg["DTOL"] = 10.0
    lib = pd.read_csv(os.path.join(HERE, cfg["library"]))
    if len(sys.argv) > 1:
        lib = lib.head(int(sys.argv[1]))

    os.makedirs(DECKS, exist_ok=True)
    manifest, skipped, inputs = [], [], []
    for _, row in lib.iterrows():
        cid = int(row["id"])
        energies = dense_energies(row["B_cm1"])
        # size RMIN against the REAL grid max (~208 cm^-1), not 64*B: feed
        # choose_rmin a one-shot energy_factor that reproduces E_max.
        emax = float(max(energies))
        rcfg = dict(cfg, energy_factors=(emax / row["B_cm1"],))
        rmin_bohr, err = mr.choose_rmin(row, rcfg)
        if err:
            skipped.append({"id": cid, "reason": err})
            continue
        deck = build_kt_deck(row, cfg, rmin_bohr, energies)
        inp = f"decks/pot_{cid:05d}.inp"
        out = f"decks/pot_{cid:05d}.out"
        with open(os.path.join(PROD, inp), "w") as fh:
            fh.write(deck)
        jtotu = mr.estimate_jtotu(max(energies), row[cfg["mu_col"]],
                                  cfg["RMAX_bohr"] * mr.BOHR2ANG)
        manifest.append({"id": cid, "rmin_bohr": round(rmin_bohr, 4),
                         "B_cm1": row["B_cm1"], "mu": row[cfg["mu_col"]],
                         "jtotu": jtotu, "n_energy": len(energies),
                         "inp": inp, "out": out})
        inputs.append(inp)

    pd.DataFrame(manifest).to_csv(os.path.join(PROD, "manifest.csv"), index=False)
    pd.DataFrame(skipped).to_csv(os.path.join(PROD, "skipped.csv"), index=False)
    with open(os.path.join(PROD, "inputs.list"), "w") as fh:
        fh.write("\n".join(inputs) + "\n")

    n = len(inputs)
    sbatch = f"""#!/bin/bash -l
#SBATCH --job-name=dv2xsec
#SBATCH --account={ACCOUNT}
#SBATCH --partition={PARTITION}
#SBATCH --array=1-{n}%{MAX_CONCURRENT}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=logs/slurm-%A_%a.log
set -e
module load imkl/2023.1.0 intel-compilers/2023.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
INP=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" inputs.list)
OUT="${{INP%.inp}}.out"
# resume: skip if already finished
if grep -q "This run used" "$OUT" 2>/dev/null; then
  echo "$INP already done"; exit 0
fi
t0=$SECONDS
{BINARY} < "$INP" > "$OUT" 2>"${{INP%.inp}}.err"
echo "$INP done in $((SECONDS-t0))s exit ${{PIPESTATUS[0]}}"
"""
    os.makedirs(os.path.join(PROD, "logs"), exist_ok=True)
    with open(os.path.join(PROD, "run_campaign.sbatch"), "w") as fh:
        fh.write(sbatch)

    print(f"runnable decks : {n}")
    print(f"skipped        : {len(skipped)} (see production/skipped.csv)")
    if skipped:
        rs = pd.DataFrame(skipped)["reason"].str.split(":").str[0].value_counts()
        for r, c in rs.items():
            print(f"    {r:32s} {c}")
    print(f"array          : 1-{n}%{MAX_CONCURRENT}")
    print(f"submit         : (cd {PROD} && sbatch run_campaign.sbatch)")


if __name__ == "__main__":
    main()
