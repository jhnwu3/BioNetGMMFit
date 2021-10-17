#!/bin/bash
#SBATCH --time=0-80:10:00 
#SBATCH --job-name=linearDebug
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/linearDebug%j.txt
#SBATCH --cpus-per-task=1
# load all modules, build terminal code, move all outputs into output folders.

./linear