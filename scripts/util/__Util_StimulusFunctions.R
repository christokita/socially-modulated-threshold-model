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
  rownames(stim) <- paste0("t", 0:gens)
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
  stim_matrix[time_step + 1, ] <- new_values #inserts stimulus into current time step (stimMat starts with row 1 = time step 0)
  return(stim_matrix)
}
