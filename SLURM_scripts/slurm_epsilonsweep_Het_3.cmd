#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 47:59:00
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=ctokita@princeton.edu

cd
Rscript scripts/3_para_sweep/3b_EpsilonParaSweep_Het_3.R