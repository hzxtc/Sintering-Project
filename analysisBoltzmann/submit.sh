#!/bin/bash

#$ -l h_data=1G,h_rt=01:00:00
#$ -cwd
#$ -o LOG.$JOB_ID
#$ -j y
#$ -M ppoths@mail
#$ -m bea
#$ -q pod_ib40.q

source /u/local/Modules/default/init/modules.sh
module load intel/16.0.2
module load intelmpi/5.0.0
module load python/3.7.2
#VASPHOME=/u/project/ana/hanguo/VASP544/vasp.5.4.4/bin/

python3 sintering.py

mv metropolis metropolis1

python3 sintering.py

mv metropolis metropolis2

python3 sintering.py

mv metropolis metropolis3

python3 sintering.py

mv metropolis metropolis4

python3 sintering.py

mv metropolis metropolis5

#date
#mkdir -p $SCRATCH/JOB.$JOB_ID
#cp -r * $SCRATCH/JOB.$JOB_ID
#echo 'SCRATCH: ' $SCRATCH/JOB.$JOB_ID
#XPWD=$PWD
#cd $SCRATCH/JOB.$JOB_ID
#mpirun -n 16 $VASPHOME/vasp_std > OUT.$JOB_ID
#cp -r * $XPWD
#date
