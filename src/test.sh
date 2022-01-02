#!/bin/bash
#SBATCH --time=0-50:10:00 
#SBATCH --job-name=parallelPSO
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=test%j.txt
#SBATCH --cpus-per-task=30

# loads necessary modules on NCH cluster, use source load.sh
git pull
module load Eigen GCC/9.3.0 OpenMPI/4.0.3 Boost/1.72.0 gnuplot
make

# load all modules, build terminal code, move all outputs into output folders.
chmod u+x sig
./sig
