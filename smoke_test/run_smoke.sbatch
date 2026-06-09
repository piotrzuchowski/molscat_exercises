#!/bin/bash -l
#SBATCH --job-name=molsmoke
#SBATCH --account=hpc-pzuch007-1692967268
#SBATCH --partition=lem-cpu-short
#SBATCH --array=1-9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=slurm-%A_%a.log

# One molscat-basic run per array task (one deck from inputs.list).
set -e
module load imkl/2023.1.0 intel-compilers/2023.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR"
BIN=/home/pzuch007/molscat_exercises/source_code/molscat-basic

INP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" inputs.list)
OUT="${INP%.input}.out"

echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: $INP -> $OUT on $(hostname)"
t0=$SECONDS
"$BIN" < "$INP" > "$OUT" 2> "${INP%.input}.err"
echo "[$(date)] done $INP in $((SECONDS-t0))s, exit ${PIPESTATUS[0]}"
