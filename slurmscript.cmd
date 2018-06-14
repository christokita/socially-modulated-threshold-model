#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 6:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=ctokita@princeton.edu

cd
Rscript scripts/2_Cluster_SocialThreshModel.R
