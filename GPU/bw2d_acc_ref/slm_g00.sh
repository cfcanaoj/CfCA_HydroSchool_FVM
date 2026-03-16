#! /bin/bash
#SBATCH --partition=workshop-1gpu
#SBATCH --gres=gpu:1
#SBATCH -o ./out%j.log
#SBATCH -e ./err%j.log


# regular partition
##SBATCH --partition=ga80-1gpu

# usage sbatch slm_g00.sh
# other useful commands
# sinfo
# squeue

module purge
module load nvhpc/25.7

./bw.x
