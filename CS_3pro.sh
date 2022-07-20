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
./CyGMM -m example/3_prot_linear_sim/3pro.bngl -c example/3_prot_linear_sim/configs/Config3pro_1.csv -x example/3_prot_linear_sim/X -y example/3_prot_linear_sim/Y -r example/3_prot_linear_sim/true_rates.csv -t example/3_prot_linear_sim/time_steps.csv
./CyGMM -m example/3_prot_linear_sim/3pro.bngl -c example/3_prot_linear_sim/configs/Config3pro_2.csv -x example/3_prot_linear_sim/X -y example/3_prot_linear_sim/Y -r example/3_prot_linear_sim/true_rates.csv -t example/3_prot_linear_sim/time_steps.csv
./CyGMM -m example/3_prot_linear_sim/3pro.bngl -c example/3_prot_linear_sim/configs/Config3pro_3.csv -x example/3_prot_linear_sim/X -y example/3_prot_linear_sim/Y -r example/3_prot_linear_sim/true_rates.csv -t example/3_prot_linear_sim/time_steps.csv
./CyGMM -m example/3_prot_linear_sim/3pro.bngl -c example/3_prot_linear_sim/configs/Config3pro_4.csv -x example/3_prot_linear_sim/X -y example/3_prot_linear_sim/Y -r example/3_prot_linear_sim/true_rates.csv -t example/3_prot_linear_sim/time_steps.csv
./CyGMM -m example/3_prot_linear_sim/3pro.bngl -c example/3_prot_linear_sim/configs/Config3pro_5.csv -x example/3_prot_linear_sim/X -y example/3_prot_linear_sim/Y -r example/3_prot_linear_sim/true_rates.csv -t example/3_prot_linear_sim/time_steps.csv

