#!/bin/bash -l
#SBATCH -A snic2021-22-850
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH -J host_job
#SBATCH --mail-type=ALL
#SBATCH --mail-user lotta.rydin@telia.com

# Load modules


# Your commands
python3 ../../Applied_Bioinformatics_OGT/final_data/collect_data_final.py