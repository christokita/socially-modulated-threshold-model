##################################################
#
# Global Stimulus Functions 
#
##################################################


####################
# Seed Stimuls
####################
seed_stimuls <- function(intitial_stim, gens) {
  # Calculate number of blank spots to make
  rep_length <- (length(intitial_stim)) * gens #intiial row does not count as gen
  # Build matrix
  stim <- matrix(data = c(intitial_stim, rep(NA, rep_length)),
                 byrow = TRUE, 
                 nrow = (gens + 1))
  # Fix Names
  colnames(stim) <- paste0(rep("s", length(intitial_stim)), 1:length(intitial_stim))
  rownames(stim) <- paste0("Gen", 0:gens)
  # Return
  return(stim)
}

####################
# Stimulus Level
####################
# Frequency depende
update_stim <- function(stim_matrix, deltas, alpha, state_matrix, time_step) {
  # Get preliminary info
  n <- nrow(state_matrix)
  stim_values <- stim_matrix[time_step, ]
  active_count <- colSums(state_matrix)
  # Calculate
  new_values <- stim_values + deltas - (alpha * (active_count/ n))
  new_values[new_values < 0] <- 0
  # Insert
  stim_matrix[time_step + 1, ] <- new_values
  return(stim_matrix)
}

# Density dependent (per capita)
globalStimUpdate_PerCap <- function(stimulus, delta, alpha, Ni, n, m, quitP) {
  # Calculate
  s <- stimulus + (alpha * (delta ) * (n / m)) - (Ni * alpha)
  # s <- stimulus + (alpha * (delta ) * (1 / (1 + quitP)) * (n / m)) - (Ni * alpha)
  # If negative, make zero
  if(s < 0.0001) {
    s <- 0
  }
  return(s)
}

####################
# Update Exhuastion
####################
updateExhaustStim <- function(ExhaustStim, TaskExhaustVect, ExhaustRate, StateMat) {
  # Zero out stim for those that are inactive
  for (i in 1:nrow(StateMat)) {
    if (sum(StateMat[i, ]) == 0) {
      ExhaustStim[i] <- 0
    }
  }
  # Get relative rate increase for task being performed
  TaskRate <- StateMat %*% TaskExhaustVect
  # Get rate increase for time step
  StimIncrease <- TaskRate %*% ExhaustRate
  # Update Stim
  NewExhaust <- ExhaustStim + StimIncrease
  # Return 
  return(NewExhaust)
}