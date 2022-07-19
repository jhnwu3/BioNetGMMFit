#!/bin/bash
#SBATCH --time=0-50:10:00 
#SBATCH --job-name=6protein
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/CS_6pro%j.txt
#SBATCH --cpus-per-task=30
# cp example/6_prot_nonlinear_slim/Y/Yt1slim.csv data/Y
# cp example/6_prot_nonlinear_slim/Y/Yt2slim.csv data/Y
# cp example/6_prot_nonlinear_slim/Y/Yt3slim.csv data/Y
# cp example/6_prot_nonlinear_slim/Y/Yt4slim.csv data/Y
./CyGMM -m 6pro.bngl -c configs/Config6pro_1.csv -x example/6_pro_nonlin_sim_slim/X/X.csv -y example/6_pro_nonlin_sim_slim/Y/Y.csv -r true_rates6.csv -t time_steps6.csv
./CyGMM -m 6pro.bngl -c configs/Config6pro_2.csv -x example/6_pro_nonlin_sim_slim/X/X.csv -y example/6_pro_nonlin_sim_slim/Y/Y.csv -r true_rates6.csv -t time_steps6.csv
./CyGMM -m 6pro.bngl -c configs/Config6pro_3.csv -x example/6_pro_nonlin_sim_slim/X/X.csv -y example/6_pro_nonlin_sim_slim/Y/Y.csv -r true_rates6.csv -t time_steps6.csv
./CyGMM -m 6pro.bngl -c configs/Config6pro_4.csv -x example/6_pro_nonlin_sim_slim/X/X.csv -y example/6_pro_nonlin_sim_slim/Y/Y.csv -r true_rates6.csv -t time_steps6.csv
./CyGMM -m 6pro.bngl -c configs/Config6pro_5.csv -x example/6_pro_nonlin_sim_slim/X/X.csv -y example/6_pro_nonlin_sim_slim/Y/Y.csv -r true_rates6.csv -t time_steps6.csv


