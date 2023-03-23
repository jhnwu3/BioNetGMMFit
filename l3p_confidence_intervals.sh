#!/bin/bash
#SBATCH --time=0-99:10:00 
#SBATCH --job-name=3proteinlinear
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/3pro%j.txt
#SBATCH --cpus-per-task=30

./BNGMM -m example/l3p_100_sim/model.bngl -x example/l3p_100_sim/X -y example/l3p_100_sim/Y -t example/l3p_100_sim/time_steps.csv -r example/l3p_100_sim/true_rates.csv -c example/l3p_100_sim/Config.csv
# ./BNGMM -m example/3_prot_linear_sim/model.bngl -x example/3_prot_linear_sim/X -y example/3_prot_linear_sim/Y -t example/3_prot_linear_sim/time_steps.csv -r example/3_prot_linear_sim/true_rates.csv -c example/3_prot_linear_sim/Config.csv