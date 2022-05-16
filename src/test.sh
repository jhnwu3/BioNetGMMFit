#!/bin/bash
#SBATCH --time=0-30:10:00 
#SBATCH --job-name=CyGMM
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=CyGMM%j.txt
#SBATCH --cpus-per-task=30

# please just load this in. :l
git pull
module load Eigen GCC/9.3.0 OpenMPI/4.0.3 Boost/1.72.0 gnuplot
make

# load all modules, build terminal code, move all outputs into output folders.
chmod u+x sig
./sig
