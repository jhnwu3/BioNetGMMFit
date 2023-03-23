#!/bin/bash
#SBATCH --time=0-50:10:00 
#SBATCH --job-name=4proteinlinear
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/4pro%j.txt
#SBATCH --cpus-per-task=30

# rm -r data/X/*
# rm -r data/Y/*
# cp example/4_prot_CD3_CD8_CD28/1min_2min/X/t1m_processed.csv data/X
# cp example/4_prot_CD3_CD8_CD28/1min_2min/Y/t2m_processed.csv data/Y
./CyGMM -m example/6_prot_nonlinear_int/6pro.bngl -c example/6_prot_nonlinear_int/Config.csv -x example/6_prot_nonlinear_int/X/ -y example/6_prot_nonlinear_int/Y/ -t example/6_prot_nonlinear_int/time_steps.csv -r example/6_prot_nonlinear_int/true_rates.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_bidir.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_ring.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_rand.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_reducedK5.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g