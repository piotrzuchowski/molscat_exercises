#!/bin/bash

name1=$1.inp
name2=$1.out
name3=$1.err
nnod=1

NPROC=8
walltime='2:00:00'
nm=$1


CURDIR=`pwd`


cat << EOF |qsub -r n
#!/bin/bash
#PBS -q main    
#PBS -l sel$ct=$nnod:ncpus=1:mpiprocs=1:mem=4000MB
#PBS -m be
#PBS -N mst_${nm:0:10}
#PBS -l $walltime

cd $CURDIR

source /usr/local/Modules/3.2.7/init/bash
module add intel
module add imkl
/home/pzuch007/NOpH2_calc/molscat/x86_64/m_Li.2S+Er.3H_wf_whf_cpld  < $name1 > $name2 


EOF



