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
./BNGMM -m 4proV2.bngl -c Config4pro.csv -x example/4_prot_CD3_CD8_CD28/1min_2min/X/ -y example/4_prot_CD3_CD8_CD28/1min_2min/Y/ -t time_steps4.csv -o frontend/output/CD8/
# ./CyGMM -m 4pro_bidir.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_ring.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_rand.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_reducedK5.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g