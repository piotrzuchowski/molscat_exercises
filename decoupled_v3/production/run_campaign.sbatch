#!/bin/bash -l
#SBATCH --job-name=dv3xsec
#SBATCH --account=hpc-pzuch007-1692967268
#SBATCH --partition=lem-cpu-short
#SBATCH --array=1-180
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
IDX=$(( ${OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
INP=$(sed -n "${IDX}p" inputs.list)
[ -z "$INP" ] && { echo "no deck at line $IDX"; exit 0; }
OUT="${INP%.inp}.out"
if grep -q "This run used" "$OUT" 2>/dev/null; then echo "$INP done"; exit 0; fi
t0=$SECONDS
/home/pzuch007/molscat_exercises/source_code/molscat-basic < "$INP" > "$OUT" 2>"${INP%.inp}.err"
echo "$INP done in $((SECONDS-t0))s exit ${PIPESTATUS[0]}"
