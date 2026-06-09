#!/bin/bash -l
#SBATCH --job-name=opacity
#SBATCH --account=hpc-pzuch007-1692967268
#SBATCH --partition=lem-cpu-short
#SBATCH --array=1-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=slurm-%A_%a.log
set -e
module load imkl/2023.1.0 intel-compilers/2023.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
INP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" inputs.list)
OUT="${INP%.input}.out"
t0=$SECONDS
/home/pzuch007/molscat_exercises/source_code/molscat-basic < "$INP" > "$OUT" 2>"${INP%.input}.err"
echo "$INP done in $((SECONDS-t0))s exit ${PIPESTATUS[0]}"
