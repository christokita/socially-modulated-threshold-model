#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=ctokita@princeton.edu

cd SocialInteractionModel/
Rscript scripts/supp_analysis/Check_Simulation_Length.R