################################################################################
#
# Social interaction model: set for running on Della cluster
#
################################################################################

rm(list = ls())

####################
# Source necessary scripts/libraries
####################
source("scripts/util/__Util__MASTER.R")
library(parallel)
library(snowfall)

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- seq(5, 100, 5) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 50000 #number of generations to run simulation 
reps           <- 100 #number of replications per simulation (for ensemble)
chunk_size     <- 5 #number of simulations sent to single core 
corrStep       <- 200 #number of time steps for calculation of correlation 

# Threshold Parameters
ThreshM        <- rep(50, m) #population threshold means 
ThreshSD       <- ThreshM * 0 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.8, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
p              <- 1 #baseline probablity of initiating an interaction per time step
epsilon        <- 0 #relative weighting of social interactions for adjusting thresholds
beta           <- 1.1 #probability of interacting with individual in same state relative to others


####################
# Prep for Parallelization
####################
# Create directory for depositing data
storage_path <- "/scratch/gpfs/ctokita/"
dir_name <- paste0("Sigma", (ThreshSD/ThreshM)[1], "-Epsilon", epsilon, "-Beta", beta, "_RankCorr")
full_path <- paste0(storage_path, dir_name)
dir.create(full_path)
sub_dirs <- c("TaskDist", "Entropy", "Thresh", "Graphs", "RankCorr")
for (sub_dir in sub_dirs) {
  dir.create(paste0(full_path, "/", sub_dir), showWarnings = FALSE)
}

# Break up parameter replications into smaller batches\
chunk_run  <- 1:(reps / chunk_size)
run_in_parallel <- expand.grid(n = Ns, run = chunk_run)
run_in_parallel <- run_in_parallel %>% 
  arrange(n)

# Prepare for parallel
no_cores <- detectCores()
sfInit(parallel = TRUE, cpus = no_cores)
sfExportAll()
sfLibrary(dplyr)
sfLibrary(reshape2)
sfLibrary(igraph)
sfLibrary(ggplot2)
sfLibrary(msm)
sfLibrary(gtools)
sfLibrary(snowfall)
sfLibrary(tidyr)
sfLibrary(stringr)
sfClusterSetupRNGstream(seed = 323)

####################
# Run ensemble simulation
####################
# Loop through group size (and chucnks)
parallel_simulations <- sfLapply(1:nrow(run_in_parallel), function(k) {
  # Set group size 
  n <- run_in_parallel[k, 1]
  chunk <- run_in_parallel[k, 2]
  # Prep lists for collection of simulation outputs from this group size
  ens_taskDist    <- list()
  ens_entropy     <- list()
  ens_thresh      <- list()
  ens_graphs      <- list()
  ens_taskCorr  <- list()
  # Run Simulations
  for (sim in 1:chunk_size) {
    ####################
    # Seed structures and intial matrices
    ####################
    # Set initial probability matrix (P_g)
    P_g <- matrix(data = rep(0, n * m), ncol = m)
    # Seed task (external) stimuli
    stimMat <- seed_stimuls(intitial_stim = InitialStim, 
                            gens = gens)
    # Seed internal thresholds
    threshMat <- seed_thresholds(n = n, 
                                 m = m, 
                                 threshold_means = ThreshM,
                                 threshold_sds = ThreshSD)
    # Start task performance
    X_g <- matrix(data = rep(0, length(P_g)), ncol = ncol(P_g))
    # Create cumulative task performance matrix
    X_tot <- X_g
    # Prep correlation step matrix
    X_prev <- matrix(data = rep(0, n * m), ncol = m)
    X_prevTot <- matrix(data = rep(0, n * m), ncol = m)
    taskCorr <- list()
    taskStep <- list()
    # Create cumulative adjacency matrix
    g_tot <-  matrix(data = rep(0, n * n), ncol = n)
    colnames(g_tot) <- paste0("v-", 1:n)
    rownames(g_tot) <- paste0("v-", 1:n)

    ####################
    # Simulate individual run
    ####################
    # Run simulation
    for (t in 1:gens) { 
      # Current timestep is actually t+1 in this formulation, because first row is timestep 0
      # Update stimuli
      stimMat <- update_stim(stim_matrix = stimMat, 
                             deltas = deltas, 
                             alpha = alpha, 
                             state_matrix = X_g, 
                             time_step = t)
      # Calculate task demand based on global stimuli
      P_g <- calc_determ_thresh(time_step        = t + 1, # first row is generation 0
                                threshold_matrix = threshMat, 
                                stimulus_matrix  = stimMat)
      # Update task performance
      X_g <- update_task_performance(task_probs   = P_g,
                                     state_matrix = X_g,
                                     quit_prob    = quitP)
      # Update social network (previously this was before probability/task update)
      g_adj <- temporalNetwork(X_sub_g = X_g,
                               prob_interact = p,
                               bias = beta)
      g_tot <- g_tot + g_adj
      # Adjust thresholds
      threshMat <- adjust_thresholds_social_capped(social_network = g_adj,
                                                   threshold_matrix = threshMat,
                                                   state_matrix = X_g,
                                                   epsilon = epsilon,
                                                   threshold_max = 2 * ThreshM[1])
      # Update total task performance profile
      X_tot <- X_tot + X_g
      # Create time step for correlation
      if (t %% corrStep == 0) {
        # Get tasks performance in correlation step
        X_step <- X_tot - X_prevTot
        # Add to ensemble list of task steps
        taskStep[[t / corrStep]] <- X_step
        # Calculate rank correlation if it is not the first step
        if(sum(X_prev) != 0) {
          # Normalize
          stepNorm <- X_step / rowSums(X_step)
          prevNorm <- X_prev / rowSums(X_prev)
          # Calculate ranks
          step_ranks <- calculateTaskRank(TaskStepMat = X_step)
          prev_ranks <- calculateTaskRank(TaskStepMat = X_prev)
          # Calculate Correlation
          rankCorr <- cor(prev_ranks, step_ranks, method = "spearman")
          # Put in list
          taskCorr[[(t / corrStep) - 1]] <- diag(rankCorr)
          names(taskCorr)[(t / corrStep) - 1] <- paste0("Gen", t)
        }
        # Update previous step total matrix
        X_prevTot <- X_tot
        # Update previous step total matrix
        X_prev <- X_step
      }
    }
    
    ####################
    # Post run calculations
    ####################
    # Calculate Entropy
    entropy <- mutualEntropy(TotalStateMat = X_tot)
    entropy <- label_parallel_runs(matrix = entropy, n = n, simulation = sim, chunk = chunk)
    # Calculate total task distribution
    totalTaskDist <- X_tot / gens
    totalTaskDist <- label_parallel_runs(matrix = totalTaskDist, n = n, simulation = sim, chunk = chunk)
    # Create thresh table
    threshMat <- label_parallel_runs(matrix = threshMat, n = n, simulation = sim, chunk = chunk)
    # Create rank correlation table
    taskCorr$replicate <- sim
    # Add total task distributions, entropy values, and graphs to lists
    ens_taskDist[[sim]]    <- totalTaskDist
    ens_entropy[[sim]]     <- entropy
    ens_thresh[[sim]]      <- threshMat 
    ens_graphs[[sim]]      <- g_tot / gens
    ens_taskCorr[[sim]]  <- taskCorr
  }
  # Bind task correlation data
  # Calculate mean correlation for each n
  runCorrs <- lapply(ens_taskCorr, function(x) {
    # Unlist
    runs <- do.call("rbind", x)
    replicate <- runs[nrow(runs), ]
    replicate <- unique(replicate)
    runs <- runs[-nrow(runs), ]
    # Calculate mean
    runMean <- matrix(data = rep(NA, m), ncol =  m)
    for (column in 1:m) {
      runMean[ , column] <- mean(runs[ , column], na.rm = TRUE)
    }
    runMean <- cbind(runMean, replicate)
    colnames(runMean) <- c(paste0("Task", 1:m), "replicate")
    return(runMean)
  })
  runCorrs <- as.data.frame(do.call("rbind", runCorrs))
  runCorrs$n <- n
  runCorrs$chunk <- chunk
  
  # Bind and write
  save_parallel_data(data = ens_taskDist, 
                     path = full_path, 
                     sub_directory = "TaskDist",
                     n = n, 
                     chunk = chunk)
  save_parallel_data(data = ens_entropy, 
                     path = full_path, 
                     sub_directory = "Entropy",
                     n = n, 
                     chunk = chunk)
  save_parallel_data(data = ens_thresh, 
                     path = full_path, 
                     sub_directory = "Thresh",
                     n = n, 
                     chunk = chunk)
  save_parallel_data(data = runCorrs, 
                     path = full_path, 
                     sub_directory = "RankCorr",
                     n = n, 
                     chunk = chunk)
  save(ens_graphs, 
       file = paste0(full_path,
                     "/Graphs/", 
                     str_pad(string = n, width = 3, pad = "0"), 
                     "-", 
                     str_pad(string = chunk, width = 2, pad = "0"), 
                     ".Rdata"))
  # Return all_clear
  rm(ens_taskDist, ens_entropy, ens_thresh, ens_graphs)
  return(paste0("DONE: n = ", n, 
                ", Sims ", ((chunk-1) * chunk_size)+1, "-", chunk * chunk_size))
  sys.sleep(1)
})

sfStop()

