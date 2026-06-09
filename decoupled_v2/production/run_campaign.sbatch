#!/bin/bash -l
#SBATCH --job-name=dv2xsec
#SBATCH --account=hpc-pzuch007-1692967268
#SBATCH --partition=lem-cpu-short
#SBATCH --array=1-946%180
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=logs/slurm-%A_%a.log
set -e
module load imkl/2023.1.0 intel-compilers/2023.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
INP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" inputs.list)
OUT="${INP%.inp}.out"
# resume: skip if already finished
if grep -q "This run used" "$OUT" 2>/dev/null; then
  echo "$INP already done"; exit 0
fi
t0=$SECONDS
/home/pzuch007/molscat_exercises/source_code/molscat-basic < "$INP" > "$OUT" 2>"${INP%.inp}.err"
echo "$INP done in $((SECONDS-t0))s exit ${PIPESTATUS[0]}"
