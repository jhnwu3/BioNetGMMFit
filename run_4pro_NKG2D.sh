#!/bin/bash
#SBATCH --time=0-50:10:00 
#SBATCH --job-name=4proteinlinear
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/4pro%j.txt
#SBATCH --cpus-per-task=30

./CyGMM -m 4pro.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
./CyGMM -m 4pro_bidir.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
./CyGMM -m 4pro_ring.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
./CyGMM -m 4pro_rand.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
./CyGMM -m 4pro_reducedK5.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g