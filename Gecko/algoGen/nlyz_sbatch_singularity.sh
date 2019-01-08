#!/bin/bash
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 1-10:00              # Runtime in D-HH:MM
#SBATCH --mem=60000
#SBATCH -o plot_analysis_%j.out      # File to which STDOUT will be written
#SBATCH -e plot_analysis_%j.err      # File to which STDERR will be written



module load system/Python-3.6.3

echo $1 , $2, $3 ,$4 ,$5

singularity exec --pwd /GECKO/ GECKO python3 plotter_for_eachhistorylog.py $1 $2 $3 $4 $5 $6 >$1/nlyzlog  2>&1

