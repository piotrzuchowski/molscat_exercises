#!/usr/bin/env python3
"""Build production MOLSCAT decks for the full library, batched in groups of 50.

- Energy grid per system from energy_kinteg_helper.energy_grid(B, jmax)
  (40 log pts 0.05-50 + 12 lin to E_max=6*kB*T_hi, plus rotor-threshold brackets).
- jmax = jmax_calc.set_jmax(B, E_max, cap=JMAX_CAP)   (campaign cap = 15)
- Numerics (validated): DR=0.04, RMAX=40 Ang, DTOL=10, JSTEP=1, IRMSET=0.
- Partners: He and H2  -> 2 decks/potential.
- Layout: production/batchNN/  (50 potentials x 2 partners = 100 decks each),
  each with inputs.list + run_batch.sbatch; plus submit_all.sh and manifest.csv.

Reuses the deck-construction core in make_decks.py (smoke_test/).
"""

import os
import sys
import pandas as pd

REPO = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "smoke_test"))
import make_decks as md            # noqa: E402
import energy_kinteg_helper as ek  # noqa: E402
import jmax_calc                   # noqa: E402

# ---- production configuration (override make_decks globals) --------------- #
md.DR = 0.04
md.RMAX_FIXED = 40.0
md.JSTEP = 1
md.DTOL = 10.0

JMAX_CAP = 15
BATCH = 50
T_HI = 50.0
PARTNERS = {"He": "mu_He", "H2": "mu_H2"}
PROD = os.path.join(REPO, "production")
ACCOUNT = "hpc-pzuch007-1692967268"
PARTITION = "lem-cpu-short"
WALLTIME = "12:00:00"
MEM = "4G"
BIN = os.path.join(REPO, "source_code", "molscat-basic")

E_MAX = 6.0 * ek.KB_CM * T_HI      # 208.5 cm^-1


def avoid_thresholds(egrid, B, jmax, rel=2e-3):
    """Move any energy sitting exactly on a rotor threshold B*j(j+1) just below it.

    MOLSCAT's RBESSK aborts ('RBESSK 0') for a channel exactly at threshold
    (kappa=0). energy_kinteg_helper.energy_grid adds points at 1.0*threshold;
    shift those to threshold*(1-rel) (channel just closed, kappa>0). The helper's
    0.98x/1.02x brackets already sample around each opening.
    """
    import numpy as np
    E = np.asarray(egrid, float).copy()
    for j in range(1, jmax + 1):
        thr = B * j * (j + 1)
        mask = np.abs(E - thr) <= 1e-6 * max(thr, 1.0)
        E[mask] = thr * (1.0 - rel)
    return np.unique(E)


def sbatch_text(ndecks, batch_id):
    return f"""#!/bin/bash -l
#SBATCH --job-name=mol_b{batch_id:02d}
#SBATCH --account={ACCOUNT}
#SBATCH --partition={PARTITION}
#SBATCH --array=1-{ndecks}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time={WALLTIME}
#SBATCH --mem={MEM}
#SBATCH --output=slurm-%A_%a.log

set -e
module load imkl/2023.1.0 intel-compilers/2023.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"

INP=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" inputs.list)
OUT="${{INP%.input}}.out"
if grep -q "This run used" "$OUT" 2>/dev/null; then
  echo "[$(date)] skip $INP (already complete)"; exit 0
fi
echo "[$(date)] $INP -> $OUT on $(hostname)"
t0=$SECONDS
{BIN} < "$INP" > "$OUT" 2> "${{INP%.input}}.err"
echo "[$(date)] done $INP in $((SECONDS-t0))s exit ${{PIPESTATUS[0]}}"
"""


def main():
    lib = pd.read_csv(os.path.join(REPO, "potentials_library.csv"))
    os.makedirs(PROD, exist_ok=True)
    ids = list(lib["id"].astype(int))
    manifest, submit_lines = [], []

    for bi, start in enumerate(range(0, len(ids), BATCH)):
        chunk = ids[start:start + BATCH]
        bdir = os.path.join(PROD, f"batch{bi:02d}")
        os.makedirs(bdir, exist_ok=True)
        inputs = []
        for pid in chunk:
            p = lib.loc[lib["id"] == pid].iloc[0]
            jmax = jmax_calc.set_jmax(p["B_cm1"], E_MAX, cap=JMAX_CAP)
            egrid = ek.energy_grid(p["B_cm1"], jmax, T_hi=T_HI)
            egrid = avoid_thresholds(egrid, p["B_cm1"], jmax)
            for partner, mucol in PARTNERS.items():
                deck, meta = md.build_deck(p, partner, float(p[mucol]),
                                           egrid, jmax)
                deck = deck.replace("(smoke test)", "(production)")
                fn = f"pot{pid:04d}_{partner}.input"
                with open(os.path.join(bdir, fn), "w") as fh:
                    fh.write(deck)
                inputs.append(fn)
                manifest.append({**meta, "batch": bi, "nnrg": len(egrid),
                                 "input": fn, "dir": f"batch{bi:02d}"})
        with open(os.path.join(bdir, "inputs.list"), "w") as fh:
            fh.write("\n".join(inputs) + "\n")
        with open(os.path.join(bdir, "run_batch.sbatch"), "w") as fh:
            fh.write(sbatch_text(len(inputs), bi))
        submit_lines.append(f"(cd batch{bi:02d} && sbatch run_batch.sbatch)")

    with open(os.path.join(PROD, "submit_all.sh"), "w") as fh:
        fh.write("#!/bin/bash\n# submit every batch array job\n"
                 + "\n".join(submit_lines) + "\n")
    os.chmod(os.path.join(PROD, "submit_all.sh"), 0o755)
    pd.DataFrame(manifest).to_csv(os.path.join(PROD, "manifest.csv"), index=False)

    nb = (len(ids) + BATCH - 1) // BATCH
    nnrg = pd.DataFrame(manifest)["nnrg"]
    print(f"Built {len(manifest)} decks in {nb} batches under {PROD}")
    print(f"  potentials: {len(ids)}  partners: {list(PARTNERS)}  "
          f"jmax cap: {JMAX_CAP}  E_max: {E_MAX:.1f} cm^-1")
    print(f"  energies/deck: min {nnrg.min()}, median {int(nnrg.median())}, "
          f"max {nnrg.max()}")
    print(f"  numerics: DR={md.DR} RMAX={md.RMAX_FIXED} DTOL={md.DTOL} "
          f"JSTEP={md.JSTEP}")
    print(f"  submit with: cd production && ./submit_all.sh")


if __name__ == "__main__":
    main()
