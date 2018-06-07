##################################################
#
# Threshold Probability of Performance 
#
##################################################


####################
# Seed task thresholds
####################
seedThresholds <- function(n, m, ThresholdMeans = NULL, ThresholdSDs = NULL) {
  # Loop through tasks and sample thresholds from normal dist
  threshMat <- lapply(1:length(ThresholdMeans), function(i) {
    threshList <- rtnorm(n = n, 
                         mean = ThresholdMeans[i], 
                         sd = ThresholdSDs[i], 
                         lower = 0)
    return(threshList)
  })
  threshMat <- do.call("cbind", threshMat)
  # Fix names
  colnames(threshMat) <- paste0("Thresh", 1:length(ThresholdMeans))
  rownames(threshMat) <- paste0("v-", 1:n)
  # Return
  return(threshMat)
}

####################
# Deterministic Threshold
####################
calcThresholdDetermMat <- function(time_step, threshold_matrix, stimulus_matrix) {
  # select proper stimulus for this time step
  no_of_stim <- ncol(threshold_matrix)
  stim_this_step <- as.matrix(stimulus_matrix[time_step, 1:no_of_stim])
  # calculate threshold probabilities for one individual
  threshold_prob <- lapply(1:ncol(threshold_matrix), function(j) {
    # Check if stimulus exceeds threshold
    thresh_result <- threshold_matrix[ , j] > stim_this_step[j]
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
adjustThresholdsSocial <- function(SocialNetwork, ThresholdMatrix, X_sub_g, epsilon) {
  # Calculate "sum" of task states/probs of neighbors
  NeighborSums <- t(SocialNetwork) %*% X_sub_g #active neighbors by category
  totalSums <- rowSums(NeighborSums) #total active neighbors
  # Loop through individuals
  for (i in 1:nrow(X_sub_g)) {
    for (j in 1:ncol(X_sub_g)) {
      activeInd <- NeighborSums[i, j]
      adjust <- epsilon * ((totalSums[i] - activeInd)  - activeInd)
      ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
      if (ThresholdMatrix[i, j] < 0) {
        ThresholdMatrix[i, j] <- 0
      } 
    }
  }
  return(ThresholdMatrix)
}

####################
# Adjust thresholds due to social interactions--with maximum level
####################
adjust_thresholds_social_capped <- function(SocialNetwork, ThresholdMatrix, X_sub_g, epsilon, ThresholdMax) {
  # Calculate "sum" of task states/probs of neighbors
  NeighborSums <- t(SocialNetwork) %*% X_sub_g #active neighbors by category
  totalSums <- rowSums(NeighborSums) #total active neighbors
  # Loop through individuals
  for (i in 1:nrow(X_sub_g)) {
    for (j in 1:ncol(X_sub_g)) {
      activeInd <- NeighborSums[i, j]
      adjust <- epsilon * ((totalSums[i] - activeInd)  - activeInd)
      ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
      if (ThresholdMatrix[i, j] < 0) {
        ThresholdMatrix[i, j] <- 0
      } else if (ThresholdMatrix[i, j] > ThresholdMax) {
        ThresholdMatrix[i, j] <- ThresholdMax
      }
    }
  }
  return(ThresholdMatrix)
}





