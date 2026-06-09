#!/bin/bash -l
#SBATCH --job-name=opac_t
#SBATCH --account=hpc-pzuch007-1692967268
#SBATCH --partition=lem-cpu-short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --output=%x.log
set -e
module load imkl/2023.1.0 intel-compilers/2023.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
t0=$SECONDS
/home/pzuch007/molscat_exercises/source_code/molscat-basic < test457_He.input > test457_He.out 2>test457_He.err
echo "done in $((SECONDS-t0))s exit ${PIPESTATUS[0]}"
