#!/bin/bash
#SBATCH -n 14                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -p defq
#SBATCH -t 5-00:00              # Runtime in D-HH:MM
#SBATCH --mem=80000
#SBATCH -o slurmlogV2/V2_%j.out      # File to which STDOUT will be written
#SBATCH -e slurmlogV2/V2_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sylvain.barriere@igh.cnrs.fr  # Email to which notifications will be sent
module load cv-standard
#module load python/3.5.2
module load python/3.5.2-bz2

module load openmpi/psm2/2.0.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $1 $2 $3

sh prod_client_script_C++_V3.sh $1  > $2  2>&1

