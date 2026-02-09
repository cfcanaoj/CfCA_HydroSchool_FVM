#!/bin/bash
#PBS -N M1-2D
#PBS -q long
#PBS -m n
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -j oe

cd $PBS_O_WORKDIR
./m1.x
