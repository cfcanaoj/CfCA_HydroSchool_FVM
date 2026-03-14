#! /bin/bash
#SBATCH --partition=dgx-full
#SBATCH --nodes=1
#SBATCH --gpu-bind=closest
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH -o ./out%j.log
#SBATCH -e ./err%j.log

# usage sbatch slm_g00.sh
# other useful commands
# sinfo
# squeue

module purge
module load nvhpc/25.7


#export NVCOMPILER_ACC_NOTIFY=3
#export NV_ACC_DEBUG=0x800

date >& out${SLURM_JOB_ID}.log
time mpiexec -n 2 ./bw.x >> out${SLURM_JOB_ID}.log
date >> out${SLURM_JOB_ID}.log
