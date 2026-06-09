#!/bin/bash

tmpxxx=$1
nodes=1
cores=4
mem=4000
walltime=3 
queue=main
dd=`date --iso-8601=seconds`
#${input_fname:0:15} 
pname=$tmpxxx'_'$dd
dir=`pwd`

echo $pname
cat << EOF | qsub -r n
#!/bin/bash
#PBS -q $queue
#PBS -l select=$nodes:ncpus=$cores:mpiprocs=$cores:mem=${mem}MB
#PBS -l place=scatter
#PBS -m be
#PBS -N ${pname:0:15}
#PBS -l walltime=${walltime}:00:00

module add intel 
module add blas

cd "\$PBS_O_WORKDIR/"
ls
pwd
/home/pzuch007/molscat_exercises/source_code/molscat-Sr_O2  < $tmpxxx >$tmpxxx'.out'
EOF
