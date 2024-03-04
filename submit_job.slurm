#!/bin/bash
#SBATCH --job-name=te_to
#SBATCH --partition=short
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=07:59:59
#SBATCH --constraint=ib
#SBATCH --exclude=d0012,d0020,d0022,d0024,d0045,d0078,d0080,d0083,d0084,d0085,d0086,d0087,d0100,d0104,d0108,d0128,d0130,d0136

module load gcc/10.1.0
module load openmpi/4.0.5-skylake-gcc10.1

mpirun ./Run_Simulation $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14}

echo $SLURM_JOB_NODELIST
