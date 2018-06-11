################################################################################
#
# Test the specilization fit (% increase) over sigma and n-slope parameter space
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")
library(parallel)
library(snowfall)
library(rlecuyer)

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- c(2, 16) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 10000 #number of generations to run simulation 
corrStep       <- 200 #number of time steps for calculation of correlation 
reps           <- 100 #number of replications per simulation (for ensemble) !!Change!!

# Threshold Parameters
ThreshM        <- c(10, 10) #population threshold means 
sigmas         <- seq(0, 0.5, 0.01)
sigmasOld      <- c(0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.4, 0.45, 0.5)
sigmas         <- sigmas[!sigmas %in% sigmasOld] 
# sigmas         <- c(0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.4, 0.45, 0.5)
InitialStim    <- c(0, 0) #intital vector of stimuli
StimRates      <- c(0.6, 0.6) #vector of stimuli increase rates  
threshSlopes   <- seq(1, 30, 1) #exponent parameter for threshold curve shape
threshSlopesOld<- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 22, 25, 27, 30) #exponent parameter for threshold curve shape
threshSlopes   <- threshSlopes[!threshSlopes %in% threshSlopesOld]
# threshSlopes   <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 22, 25, 27, 30) #exponent parameter for threshold curve shape
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
p              <- 0 #probability of interacting with individual in other states
q              <- 1 #probability of interacting with individual in same state relative to others

# Necessary to check on progress of parallel simulations
iteration <- 100 / (length(threshSlopes) * length(sigmas) * length(Ns))
writeLines("0%", con = "output/ParameterExploration/ProgressOnCurrentSimulation.txt")


####################
# Run simulation multiple times
####################
# prep mega list
# improveSpec <- list()


# Prepare for parallel
no_cores <- detectCores() - 2
sfInit( parallel = TRUE, cpus = no_cores)
sfExportAll()
sfLibrary(dplyr)
sfLibrary(reshape2)
sfLibrary(igraph)
sfLibrary(ggplot2)
sfLibrary(msm)
sfLibrary(gtools)
sfLibrary(snowfall)
sfClusterSetupRNGstream(seed = 123)

# Loop through n slopes
# for (k in 1:length(threshSlopes)) {
improveSpec <- sfLapply(1:length(threshSlopes),function(k) {
  
  # Set threshSlope 
  threshSlope <- threshSlopes[k]
  
  # Prep meta-lists for collection of group size simulations
  sigmaDiff <- data.frame(PercIncrease = NULL, SlopeIncrease = NULL, SpecLarge = NULL, SpecSmall = NULL, c = NULL, threshSlope = NULL)
  
  # Loop through sigmas
  for (z in 1:length(sigmas)) {
    
    # Set sigma 
    sigma <- sigmas[z]
    ThreshSD       <- ThreshM * sigma #population threshold standard deviations
    
    # Prep meta-lists for collection of group size simulations
    groups_taskCorr  <- list()
    groups_runStims  <- list()
    
    # Loop through group sizes
    for (i in 1:length(Ns)) {
      # Set group size
      n <- Ns[i]
      
      # Prep lists for collection of simulation outputs
      ens_taskCorr  <- list()
      ens_stimDiff  <- list()
      
      # Run Simulations
      for (sim in 1:reps) {
        
        ####################
        # Seed structures and intial matrices
        ####################
        
        # Set initial probability matrix (P_g)
        P_g <- initiateProbMatrix(n = n, m = m)
        
        # Seed task (external) stimuli
        stimMat <- seedStimuls(InitialSVector = InitialStim, 
                               RateVector = StimRates, 
                               gens = gens)
        
        # Seed internal thresholds
        threshMat <- seedThresholds(n = n, 
                                    m = m, 
                                    ThresholdMeans = ThreshM, 
                                    ThresholdSDs = ThreshSD)
        
        # Start task performance
        X_g <- matrix(data = rep(0, length(P_g)), ncol = ncol(P_g))
        
        # Create cumulative task performance matrix
        X_tot <- X_g
        
        # Prep correlation step matrix
        X_prev <- matrix(data = rep(0, n * m), ncol = m)
        X_prevTot <- matrix(data = rep(0, n * m), ncol = m)
        taskCorr <- list()
        
        ####################
        # Simulate
        ####################
        # Run simulation
        for (t in 1:gens) {
          # Update stimuli
          for (j in 1:(ncol(stimMat)/2)) {
            # update stim
            stimMat[t + 1, j] <- globalStimUpdate(stimulus = stimMat[t, j],
                                                  delta = stimMat[t, j + m], 
                                                  alpha = alpha, 
                                                  Ni = sum(X_g[ , j]), 
                                                  n = n)
            # shift down delta (rate increases)
            stimMat[t + 1, j + m] <- stimMat[t, j + m]
          }
          # Update social network
          # g_adj <- temporalNetwork(X_sub_g = X_g,
          #                          p = p,
          #                          bias = q)
          # Calculate task demand based on global stimuli
          P_g <- calcThresholdProbMat(TimeStep = t + 1, # first row is generation 0
                                      ThresholdMatrix = threshMat, 
                                      StimulusMatrix = stimMat, 
                                      nSlope = threshSlope)
          # Update task performance
          X_g <- updateTaskPerformance(P_sub_g    = P_g,
                                       TaskMat    = X_g,
                                       QuitProb   = quitP)
          # Update total task performance profile
          X_tot <- X_tot + X_g
          # Create time step for correlation
          if (t %% corrStep == 0) {
            # Get tasks performance in correlation step
            X_step <- X_tot - X_prevTot
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
        # Calculate total difference in stim over averaged window at beginning and end
        stimMat <- as.data.frame(stimMat)
        stimDiff <- data.frame(task1start = mean(stimMat$s1[101:600]), #avg over 500 window after first 100 acclimation 
                               task2start = mean(stimMat$s2[101:600]),
                               task1end   = mean(stimMat$s1[(nrow(stimMat) - 499):nrow(stimMat)]), #avg over last 500 steps
                               task2end   = mean(stimMat$s2[(nrow(stimMat) - 499):nrow(stimMat)]))
        stimDiff <- stimDiff %>% 
          summarise(Task1Diff = (task1end - task1start) / task1start,
                    Task2Diff = (task2end - task2start) / task2start)
        
        # Add total task distributions, entropy values, and graphs to lists
        ens_taskCorr[[sim]]  <- taskCorr
        ens_stimDiff[[sim]] <- stimDiff
      }
      # Calculate mean correlation for each run
      runCorrs <- lapply(ens_taskCorr, function(x) {
        # Unlist
        runs <- do.call("rbind", x)
        runs <- as.data.frame(runs)
        names(runs) <- c("Task1", "Task2")
        # Calculate mean
        runs <- 
          runs %>% 
          summarise(Task1 = mean(Task1), Task2 = mean(Task2))
      })
      runCorrs <- do.call("rbind", runCorrs)
      runCorrs$n  <- n
      
      # Calculate mean stim difference
      runStims <- do.call("rbind", ens_stimDiff)
      runStims$n <- n
      
      # Add to list of lists
      groups_taskCorr[[i]]  <- runCorrs
      groups_runStims[[i]]  <- runStims
      
      #rest
      Sys.sleep(2)
      
      # Print simulation completed
      iterationOut <- readLines("output/ParameterExploration/ProgressOnCurrentSimulation.txt")
      iterationOut <- as.numeric(gsub("[^\\.0-9]", "", iterationOut))
      iterationOut <- iterationOut + iteration
      iterationOut <- round(iterationOut, digits = 2)
      writeLines(paste0(iterationOut, "%"), con = "output/ParameterExploration/ProgressOnCurrentSimulation.txt")
    }
    # Unlist
    taskCorrTot <- do.call("rbind", groups_taskCorr)
    taskCorrTot[is.na(taskCorrTot)] <- 0 #fix NAs
    taskCorrTot$TaskMean <- (taskCorrTot$Task1 +  taskCorrTot$Task2) / 2
    # Calculate mean task correlation
    taskCorrMeans <- data.frame(n = unique(taskCorrTot$n), SpecMean = NA)
    for (i in 1:nrow(taskCorrMeans)) {
      taskCorrMeans$SpecMean[i] <- mean(taskCorrTot$TaskMean[taskCorrTot$n == taskCorrMeans$n[i]])
    }
    # Unlist stim differences
    taskStimDiff <- do.call("rbind", groups_runStims)
    taskStimMean <- data.frame(n = unique(taskStimDiff$n), Task1Diff = NA, Task2Diff = NA)
    for (i in 1:nrow(taskCorrMeans)) {
      taskStimMean$Task1Diff[i] <- mean(taskStimDiff$Task1Diff[taskStimDiff$n == taskStimMean$n[i]])
      taskStimMean$Task2Diff[i] <- mean(taskStimDiff$Task2Diff[taskStimDiff$n == taskStimMean$n[i]])
    }
    # Calculate improvement in specialization between group size 2 and 16
    CorrMeanDiff <- data.frame(PercIncrease = NA, 
                               SlopeIncrease = NA, 
                               SpecSmall = NA, 
                               SpecLarge = NA, 
                               Task1DiffSmall = taskStimMean$Task1Diff[taskStimMean$n == 2],
                               Task2DiffSmall = taskStimMean$Task2Diff[taskStimMean$n == 2],
                               Task1DiffLarge = taskStimMean$Task1Diff[taskStimMean$n == 16],
                               Task2DiffLarge = taskStimMean$Task2Diff[taskStimMean$n == 16],
                               sigma = sigma, 
                               threshSlope = threshSlope)
    n2 <- taskCorrMeans$SpecMean[taskCorrMeans$n == 2]
    n16 <- taskCorrMeans$SpecMean[taskCorrMeans$n == 16]
    CorrMeanDiff$PercIncrease <- (n16 - n2) / abs(n2)
    CorrMeanDiff$SlopeIncrease <- (n16 - n2) / (16 - 2)
    CorrMeanDiff$SpecLarge <- n16
    CorrMeanDiff$SpecSmall <- n2
    # Bind to list within c value
    sigmaDiff <- rbind(sigmaDiff, CorrMeanDiff)
  }
  
  # add results to list
  return(sigmaDiff)
  
})

sfStop()
####################
# Prep and Plot
####################
# Unlist
improve <- do.call("rbind", improveSpec)
