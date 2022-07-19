#!/bin/bash
#SBATCH --time=0-50:10:00 
#SBATCH --job-name=4proteinlinear
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/CS_CD8%j.txt
#SBATCH --cpus-per-task=30

./CyGMM -m 4pro.bngl -c configs/ConfigCD8_1.csv -x example/4_prot_CD3_CD8_CD28/1min_2min/X/t1m_processed.csv -y example/4_prot_CD3_CD8_CD28/1min_2min/Y/t2m_processed.csv -t time_steps4.csv -o frontend/output/ --g
./CyGMM -m 4pro.bngl -c configs/ConfigCD8_2.csv -x example/4_prot_CD3_CD8_CD28/1min_2min/X/t1m_processed.csv -y example/4_prot_CD3_CD8_CD28/1min_2min/Y/t2m_processed.csv -t time_steps4.csv -o frontend/output/ --g
./CyGMM -m 4pro.bngl -c configs/ConfigCD8_3.csv -x example/4_prot_CD3_CD8_CD28/1min_2min/X/t1m_processed.csv -y example/4_prot_CD3_CD8_CD28/1min_2min/Y/t2m_processed.csv -t time_steps4.csv -o frontend/output/ --g
./CyGMM -m 4pro.bngl -c configs/ConfigCD8_4.csv -x example/4_prot_CD3_CD8_CD28/1min_2min/X/t1m_processed.csv -y example/4_prot_CD3_CD8_CD28/1min_2min/Y/t2m_processed.csv -t time_steps4.csv -o frontend/output/ --g
./CyGMM -m 4pro.bngl -c configs/ConfigCD8_5.csv -x example/4_prot_CD3_CD8_CD28/1min_2min/X/t1m_processed.csv -y example/4_prot_CD3_CD8_CD28/1min_2min/Y/t2m_processed.csv -t time_steps4.csv -o frontend/output/ --g



# ./CyGMM -m 4pro_bidir.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_ring.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_rand.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g
# ./CyGMM -m 4pro_reducedK5.bngl -c Config4pro.csv -t time_steps4.csv -o frontend/output/ --g