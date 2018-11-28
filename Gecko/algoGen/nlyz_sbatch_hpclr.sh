#!/bin/bash
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -p defq
#SBATCH -t 1-10:00              # Runtime in D-HH:MM
#SBATCH --mem=60000
#SBATCH -o slurmlogV2/nlyz_%j.out      # File to which STDOUT will be written
#SBATCH -e slurmlogV2/nlyz_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sylvain.barriere@igh.cnrs.fr  # Email to which notifications will be sent
module load cv-standard

module load python/3.5.2-bz2

echo $1 , $2, $3 ,$4 ,$5

python3 plotter_for_eachhistorylog.py $1 $2 $3 $4 $5 $6 >$1/nlyzlog  2>&1
