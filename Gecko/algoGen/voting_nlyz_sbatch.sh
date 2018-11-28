#!/bin/bash
#SBATCH -n 3                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH --mem=60000
#SBATCH -o voting_%j.out      # File to which STDOUT will be written
#SBATCH -e voting_%j.err      # File to which STDERR will be written


module load system/Python-3.6.3

echo $1 , $2


srun -n1  python3 occurencemerge.py $1 >$1occurencemerge.log  2>&1 &
srun -n1   python3 cmpfastqmine.py $1 >$1cmpfastqmine.log  2>&1 &
srun -n1  python3 clustering_winner.py $1 $2 >$1clustering_winner.log  2>&1 &
wait
