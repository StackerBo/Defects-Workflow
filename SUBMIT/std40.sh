#!/bin/bash
# Torque script for Wuzhou
#PBS -N O2
#PBS -q std40
#PBS -l nodes=4:ppn=40
#PBS -V
#PBS -o v1_out.log
#PBS -e v2_err.log
#PBS -l walltime=480:00:00

module load intel/intel2020
module load vasp/vasp544elpa
#module load qe/qe-6.5
#module load misc/atat/atatmod



cd $PBS_O_WORKDIR
nproc=`cat $PBS_NODEFILE | wc -l`
export MPIRUN=`echo "mpirun -np $nproc -genv I_MPI_FABRICS=shm:ofi" `



echo $PBS_JOBID >> "v_JOBID=[${PBS_JOBID}]"
#echo "$NODE nodes in $part, vasp544_$vers" | tee node_part_version.log



date >>time_start.log
$MPIRUN vasp_std >> v0_vasp.log
date >>time_finish.log
