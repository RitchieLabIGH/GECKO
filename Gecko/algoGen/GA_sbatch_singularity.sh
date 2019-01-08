#!/bin/bash
#SBATCH -n 10                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH --mem=90000
#SBATCH -o GECKO_%j.out      # File to which STDOUT will be written
#SBATCH -e GECKO_%j.err      # File to which STDERR will be written


module load compiler/gcc-5.3.0
module load mpi/openmpi-2.1.2
module load system/Python-3.6.3
export OMP_NUM_THREADS=$SLURM_NTASKS
echo $1 $2 $3

singularity exec -w --pwd /GECKO/ GECKO sh prod_client_script_C++_V3.sh $1  > $2  2>&1


