#!/bin/bash
#SBATCH --time=0-99:10:00 
#SBATCH --job-name=l6p
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/6linpro%j.txt
#SBATCH --cpus-per-task=30

# ./BNGMM -m example/l3p_100_sim/model.bngl -x example/l3p_100_sim/X -y example/l3p_100_sim/Y -t example/l3p_100_sim/time_steps.csv -r example/l3p_100_sim/true_rates.csv -c example/l3p_100_sim/Config.csv
# ./BNGMM -m example/3_prot_linear_sim/model.bngl -x example/3_prot_linear_sim/X -y example/3_prot_linear_sim/Y -t example/3_prot_linear_sim/time_steps.csv -r example/3_prot_linear_sim/true_rates.csv -c example/3_prot_linear_sim/Config.csv
# ./BNGMM -m example/yeast/yeast.bngl -x example/yeast/X -y example/yeast/Y -t example/yeast/time_steps.csv -r example/yeast/true_rates.csv -c example/yeast/Config.csv -hr example/yeast/heldRates.csv -f example/yeast/future_times.csv -o test/
# ./BNGMM -m example/6_pro_lin_sim/6_pro_lin.bngl -c example/6_pro_lin_sim/Config.csv -x example/6_pro_lin_sim/X -y example/6_pro_lin_sim/Y -t example/6_pro_lin_sim/ts.csv -r example/6_pro_lin_sim/tr.csv -o test/l6p5p_m/ --g
./BNGMM -m example/l6p_5p_sim/6_pro_lin.bngl -c example/l6p_5p_sim/ConfigM.csv -x example/l6p_5p_sim/X -y example/l6p_5p_sim/Y -t example/l6p_5p_sim/ts.csv -r example/l6p_5p_sim/tr.csv -o test/l6p5p_m/ --g
# ./BNGMM -m example/l6p_5p_sim/6_pro_lin.bngl -c example/l6p_5p_sim/ConfigMV.csv -x example/l6p_5p_sim/X -y example/l6p_5p_sim/Y -t example/l6p_5p_sim/ts.csv -r example/l6p_5p_sim/tr.csv -o test/l6p5p_mv/ --g
# ./BNGMM -m example/l6p_5p_sim/6_pro_lin.bngl -c example/l6p_5p_sim/ConfigMVC.csv -x example/l6p_5p_sim/X -y example/l6p_5p_sim/Y -t example/l6p_5p_sim/ts.csv -r example/l6p_5p_sim/tr.csv -o test/l6p5p_mvc/ --g