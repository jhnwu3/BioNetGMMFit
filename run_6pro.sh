#!/bin/bash
#SBATCH --time=0-50:10:00 
#SBATCH --job-name=6protein
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/6pro%j.txt
#SBATCH --cpus-per-task=30
rm -r data/X/*
rm -r data/Y/*
cp example/6_protein_nonlinear_sim/X/X.csv data/X
cp example/6_protein_nonlinear_sim/Y/Y.csv data/Y
./CyGMM -m 6pro.bngl -c Config6pro.csv -r true_rates6.csv -t time_steps6.csv
