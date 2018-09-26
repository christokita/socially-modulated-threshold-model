##################################################
#
# Threshold Probability of Performance 
#
##################################################

####################
# Seed task thresholds
####################
seed_thresholds <- function(n, m, threshold_means = NULL, threshold_sds = NULL) {
  # Loop through tasks and sample thresholds from normal dist
  thresh_mat <- lapply(1:length(threshold_means), function(i) {
    threshList <- rtnorm(n = n, 
                         mean = threshold_means[i], 
                         sd = threshold_sds[i], 
                         # lower = 0)
                         lower = -0.00001)
    return(threshList)
  })
  thresh_mat <- do.call("cbind", thresh_mat)
  # Fix names
  colnames(thresh_mat) <- paste0("Thresh", 1:length(threshold_means))
  rownames(thresh_mat) <- paste0("v-", 1:n)
  # Return
  return(thresh_mat)
}

####################
# Deterministic Threshold
####################
calc_determ_thresh <- function(time_step, threshold_matrix, stimulus_matrix) {
  # select proper stimulus for this time step
  no_of_stim <- ncol(threshold_matrix)
  stim_this_step <- as.matrix(stimulus_matrix[time_step, 1:no_of_stim])
  # calculate threshold probabilities for one individual
  threshold_prob <- lapply(1:ncol(threshold_matrix), function(j) {
    # Check if stimulus exceeds threshold
    thresh_result <-  stim_this_step[j] > threshold_matrix[ , j]
    return(as.numeric(thresh_result))
  })
  # bind and return
  threshold_prob <- do.call("cbind", threshold_prob)
  colnames(threshold_prob) <- paste0("ThreshProb", 1:ncol(threshold_prob))
  rownames(threshold_prob) <- paste0("v-", 1:nrow(threshold_prob))
  return(threshold_prob)
}

####################
# Adjust thresholds due to social interactions
####################
adjust_thresh_social <- function(social_network, threshold_matrix, state_matrix, epsilon) {
  # Calculate "sum" of task states/probs of neighbors
  active_neighbors <- t(social_network) %*% state_matrix #active neighbors by category
  total_neighbors <- rowSums(active_neighbors) #total active neighbors
  # Calculate interacting individuals NOT performing that task
  not_sums <- total_neighbors - active_neighbors
  # Calculate net threshold effect
  net_effect <- not_sums - active_neighbors
  net_effect <- net_effect * epsilon
  # Adjust thresholds and return
  threshold_matrix <- threshold_matrix + net_effect
  # Catch minimum thresholds
  threshold_matrix[threshold_matrix < 0] <- 0
  # Return
  return(threshold_matrix)
}

####################
# Adjust thresholds due to social interactions--with maximum level
####################
adjust_thresholds_social_capped <- function(social_network, threshold_matrix, state_matrix, epsilon, threshold_max) {
  # Calculate "sum" of task states/probs of neighbors
  active_neighbors <- t(social_network) %*% state_matrix #active neighbors by category
  total_neighbors <- rowSums(active_neighbors) #total active neighbors
  # Calculate interacting individuals NOT performing that task
  not_sums <- total_neighbors - active_neighbors
  # Calculate net threshold effect
  net_effect <- not_sums - active_neighbors
  net_effect <- net_effect * epsilon
  # Adjust thresholds and return
  threshold_matrix <- threshold_matrix + net_effect
  # Catch minimum and maximum thresholds
  threshold_matrix[threshold_matrix < 0] <- 0
  threshold_matrix[threshold_matrix > threshold_max] <- threshold_max
  # Return
  return(threshold_matrix)
}

####################
# Summarize threshold tracking matrix
####################
summarise_threshold_tracking <- function(tracked_threshold, n, time_steps, chunk, simulation) {
  # Get column names
  id_names <- names(tracked_threshold[[1]])
  # Unlist
  threshtime <- matrix(unlist(tracked_threshold),
                       ncol = length(tracked_threshold[[1]]),
                       byrow = TRUE,
                       dimnames = list(c(NULL), c(id_names)))
  # Reshape
  row.names(threshtime) <- NULL
  threshtime <- as.data.frame(threshtime)
  threshtime <- threshtime %>% 
    gather("Id", "Threshold") 
  threshtime$t <- rep(0:time_steps, n)
  threshtime$n <- n
  threshtime$simulation <- simulation
  threshtime$chunk <- chunk
  # Return
  return(threshtime)
}




