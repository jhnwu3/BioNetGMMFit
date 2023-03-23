#!/bin/bash
#SBATCH --time=0-50:10:00 
#SBATCH --job-name=6protein
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/6pro%j.txt
#SBATCH --cpus-per-task=30
rm -r data/X/*
rm -r data/Y/*
cp example/6_pro_nonlin_sim_slim/X/X.csv data/X
cp example/6_pro_nonlin_sim_slim/Y/Y.csv data/Y
# cp example/6_prot_nonlinear_slim/Y/Yt1slim.csv data/Y
# cp example/6_prot_nonlinear_slim/Y/Yt2slim.csv data/Y
# cp example/6_prot_nonlinear_slim/Y/Yt3slim.csv data/Y
# cp example/6_prot_nonlinear_slim/Y/Yt4slim.csv data/Y
./BNGMM -m 6pro.bngl -c Config6pro.csv -r true_rates6.csv -t time_steps6.csv
