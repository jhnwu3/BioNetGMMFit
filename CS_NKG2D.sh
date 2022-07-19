#!/bin/bash
#SBATCH --time=0-50:10:00 
#SBATCH --job-name=NKG2D
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurm_outputs/CS_NKG2D%j.txt
#SBATCH --cpus-per-task=30

./CyGMM -m NKG2D.bngl -c configs/ConfigNKG2D_1.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
./CyGMM -m NKG2D.bngl -c configs/ConfigNKG2D_2.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
./CyGMM -m NKG2D.bngl -c configs/ConfigNKG2D_3.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
./CyGMM -m NKG2D.bngl -c configs/ConfigNKG2D_4.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
./CyGMM -m NKG2D.bngl -c configs/ConfigNKG2D_5.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g




# ./CyGMM -m 4pro_bidir.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
# ./CyGMM -m 4pro_ring.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
# ./CyGMM -m 4pro_rand.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g
# ./CyGMM -m 4pro_reducedK5.bngl -c Config4pro.csv -t time_steps4nkg2d.csv -o frontend/output/ -x example/4_pro_NKG2D_dim/X -y example/4_pro_NKG2D_dim/Y --g