#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=M-large-t
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem=30G
#SBATCH -o %x_%J.out
#SBATCH -e %x_%J.out
#SBATCH --hint=nomultithread
#SBATCH --reservation=hydro

### the number of threads (can be modified)
#SBATCH --cpus-per-task=32

module switch PrgEnv-cray PrgEnv-intel
module list

export KMP_AFFINITY=compact

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} 

cd $SLURM_SUBMIT_DIR

srun -c ${SLURM_CPUS_PER_TASK} ./a.out
