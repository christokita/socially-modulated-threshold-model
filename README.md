# Socially Modulated Threshold Model
Agent-based model investigating how social interactions can modify response thresholds and result in self-organized behavioral specialization (i.e, division of labor) and social network structure.


## Overview
Scripts of the computational model and subsequent simulation data for:

>Tokita CK & Tarnita CE. (2020). Social influence and interaction bias can drive emergent behavioural specialization and modular social networks across systems. *Journal of the Royal Society Interface* 17: 20190564. [http://dx.doi.org/10.1098/rsif.2019.0564](http://dx.doi.org/10.1098/rsif.2019.0564)

This is the living repository of this code. For the public release/archive of the code used in the paper please see [release v1.2.0](https://github.com/christokita/socially-modulated-thresholds/releases/tag/v1.2.0) of this repository.

## Components of this repository

This directory has three main files:

* **scripts**: contains all scripts for simulation, analysis, and plotting of the computational model and analytical model.
* **output**: contains all subsequent derived data from simulations as well as any graphs produced from analysis of the simulation data.
* **SLURM_scripts**: contains all scripts for running R code/simulations on Princeton Della Clusters. I recently reorganized files for this release, so be careful with file names/paths in the scripts here. 

Scripts are numbered numerically for organization. The scripts folder has the following structure:
* **util** folder contains all custom-made functions used throughout the simualtion. This then, in essence, contains the model itself.
* **0** installs packages necessary for running the model in parallel on computer cluster. Run this on head node of computer cluster first.
* **1** non-parallel runs of the model.
* **2** Parallelized runs of the model. 
  * **2b** allows processing of data produced by these scripts. 
* **3** Parallelized runs of the model for large parameter sweeps across interaction bias (beta) and social influence (epsilon). 
  * **3b** can allows processing of data produced by simulations. 
  * **3c** allows plotting of this processed data.
  * **3_para_sweep** folder contains extensive scripts for breaking up parameter sweeps into smaller sweeps for faster simulation on multiple nodes.
* **Analyze_** scripts are used for analyzing various simulation data, organized by focus of the analysis.
  * **Analyze_DOL_** scripts focus on analyzing task-performance/behavior data.
  * **Analyze_Network_** scripts focus on analyzing interactions patterns, i.e., time-aggregated social networks. 
  * **Analyze_Stimulus** analyzes task-need stimuli over time.
  * **Analyze_TaskDistribution** analyzes the overall frequency of task performance in a run of the model (usually sweeping across one parameter). 
  * **Analyze_ThresholdDistribution** analyzes the distribution of final thresholds in a run of the model (usually sweeping across one parameter). 
  * **Analyze_ThreshOverTime** analyzes the time evolution of thresholds in a single run of a model (usually sweeping across group size). 
* **analytical_calcs** folder contains R and Mathematica scripts for the analytical calculations detailed in the supplemental information. 
* **long_sims** folder contains scripts for running simulations of the model on much longer time scales and analyzing effect of long simulation on model results. 
* **supp_analysis** folder contains various scripts for checking various details of the model (e.g., changing interaction probability, removing threshold limit)
* **testing_scripts** folder contains various "sanity-check" simulations to make sure parallelization does not affect random number generator in model. Also, contains a way to rename files. 
