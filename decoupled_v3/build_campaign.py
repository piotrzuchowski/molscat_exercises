#!/usr/bin/env python3
"""Build the decoupled_v3 cross-section / k(T) campaign for the cluster.

Uses molscat_run.build_deck (R^-7 lambda=1 term) + the jmax15 dense k(T) energy
grid. Writes one deck per potential, a manifest, inputs.list and a SLURM array
script. Run the chunked controller (run_all.sh) to submit within the 200 cap.

    python build_campaign.py [n_cases]      # default: all rows in the library
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
import molscat_run as mr
import energy_kinteg_helper as ek

BINARY = mr.CONFIG["binary"]
ACCOUNT = "hpc-pzuch007-1692967268"
PARTITION = "lem-cpu-short"
THROTTLE = 180
JMAX = mr.CONFIG["JMAX"]
T_HI = 50.0
PROD = os.path.join(HERE, "production")
DECKS = os.path.join(PROD, "decks")


def avoid_thresholds(egrid, B, jmax, rel=2.0e-3):
    E = np.asarray(egrid, float).copy()
    for j in range(1, jmax + 1):
        t = B * j * (j + 1)
        if t <= 0:
            continue
        E[np.abs(E - t) <= rel * t] = t * (1.0 - rel)
    return np.unique(E)


def dense_energies(B):
    return avoid_thresholds(ek.energy_grid(B, JMAX, T_hi=T_HI), B, JMAX)


def main():
    cfg = dict(mr.CONFIG)
    lib = pd.read_csv(os.path.join(HERE, cfg["library"]))
    if len(sys.argv) > 1:
        lib = lib.head(int(sys.argv[1]))
    os.makedirs(DECKS, exist_ok=True)
    os.makedirs(os.path.join(PROD, "logs"), exist_ok=True)

    manifest, skipped, inputs = [], [], []
    for _, row in lib.iterrows():
        cid = int(row["id"])
        energies = dense_energies(row["B_cm1"])
        cfg["_emax_cm1"] = float(max(energies))
        rmin_bohr, err = mr.choose_rmin(row, cfg)
        if err:
            skipped.append({"id": cid, "reason": err})
            continue
        deck = mr.build_deck(row, cfg, rmin_bohr, energies)
        inp = f"decks/pot_{cid:05d}.inp"
        with open(os.path.join(PROD, inp), "w") as fh:
            fh.write(deck)
        jtotu = mr.estimate_jtotu(cfg["_emax_cm1"], row[cfg["mu_col"]],
                                  cfg["RMAX_bohr"] * mr.BOHR2ANG)
        manifest.append({"id": cid, "rmin_bohr": round(rmin_bohr, 4),
                         "B_cm1": row["B_cm1"], "mu": row[cfg["mu_col"]],
                         "jtotu": jtotu, "n_energy": len(energies), "inp": inp,
                         "out": inp.replace(".inp", ".out")})
        inputs.append(inp)

    pd.DataFrame(manifest).to_csv(os.path.join(PROD, "manifest.csv"), index=False)
    pd.DataFrame(skipped).to_csv(os.path.join(PROD, "skipped.csv"), index=False)
    with open(os.path.join(PROD, "inputs.list"), "w") as fh:
        fh.write("\n".join(inputs) + "\n")

    n = len(inputs)
    sbatch = f"""#!/bin/bash -l
#SBATCH --job-name=dv3xsec
#SBATCH --account={ACCOUNT}
#SBATCH --partition={PARTITION}
#SBATCH --array=1-{THROTTLE}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --mem=2G
#SBATCH --output=logs/slurm-%A_%a.log
set -e
module load imkl/2023.1.0 intel-compilers/2023.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
# deck index = OFFSET + array task id (array ids <= MaxArraySize=1001)
IDX=$(( ${{OFFSET:-0}} + SLURM_ARRAY_TASK_ID ))
INP=$(sed -n "${{IDX}}p" inputs.list)
[ -z "$INP" ] && {{ echo "no deck at line $IDX"; exit 0; }}
OUT="${{INP%.inp}}.out"
if grep -q "This run used" "$OUT" 2>/dev/null; then echo "$INP done"; exit 0; fi
t0=$SECONDS
{BINARY} < "$INP" > "$OUT" 2>"${{INP%.inp}}.err"
echo "$INP done in $((SECONDS-t0))s exit ${{PIPESTATUS[0]}}"
"""
    with open(os.path.join(PROD, "run_campaign.sbatch"), "w") as fh:
        fh.write(sbatch)

    print(f"runnable decks : {n}")
    print(f"skipped        : {len(skipped)}")
    if skipped:
        for r, c in pd.DataFrame(skipped)["reason"].str.split(":").str[0]\
                .value_counts().items():
            print(f"    {r:28s} {c}")
    print(f"array 1-{n}%{THROTTLE}  ->  run via run_all.sh controller")


if __name__ == "__main__":
    main()
