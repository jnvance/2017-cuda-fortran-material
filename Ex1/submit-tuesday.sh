#!/bin/bash
#SBATCH -J NAME
#SBATCH -A TRAINING-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:30:00
#SBATCH --no-requeue
#SBATCH --partition=tesla
#SBATCH --reservation=cuda_fortran_tue

# EDIT ME ONLY IF YOU KNOW WHAT YOU ARE DOING ###############
. /etc/profile.d/modules.sh
module purge
module load default-wilkes
module unload cuda intel/cce intel/fce intel/impi intel/mkl
module load intel/mkl/11.3.3.210
module load pgi/16.10
module load openmpi/pgi/1.10.3
module load cuda/8.0
module load custom/magma/2.0.2 custom/lib-jdr/1.0
#############################################################

cd $SLURM_SUBMIT_DIR

echo "Adding with CPU"
./vadd_cpu
echo "Adding with GPU"
./vadd_gpu
