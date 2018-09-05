#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 1:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=ctokita@princeton.edu

cd
Rscript scripts/testing_scripts/FixedThresh_VariationSweep.R