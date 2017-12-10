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
calcThresholdDetermMat <- function(TimeStep, ThresholdMatrix, StimulusMatrix) {
  # select proper stimulus for this time step
  stimulusThisStep <- StimulusMatrix[TimeStep, ]
  # calculate threshold probabilities for one individual
  thresholdP <- lapply(1:nrow(ThresholdMatrix), function(i) {
    # select row for individual in threshold matrix
    indThresh <- ThresholdMatrix[i, ]
    # create task vector to be output and bound
    taskThresh <- rep(0, length(indThresh))
    # loop through each task within individual
    for (j in 1:length(taskThresh)) {
      stim <- stimulusThisStep[j]
      thresh <- indThresh[j]
      if (stim > thresh) {
        taskThresh[j] <- 1
      }
    }
    return(taskThresh)
  })
  # bind and return
  thresholdP <- do.call("rbind", thresholdP)
  colnames(thresholdP) <- paste0("ThreshProb", 1:ncol(thresholdP))
  rownames(thresholdP) <- paste0("v-", 1:nrow(thresholdP))
  return(thresholdP)
}

####################
# Deterministic Threshold with social information
####################
calcThresholdDetermSocial <- function(TimeStep, ThresholdMatrix, StimulusMatrix, epsilon, SocialInfoMatrix) {
  # select proper stimulus for this time step
  stimulusThisStep <- StimulusMatrix[TimeStep, ]
  # calculate threshold probabilities for one individual
  thresholdP <- lapply(1:nrow(ThresholdMatrix), function(i) {
    # select row for individual in threshold matrix
    indThresh <- ThresholdMatrix[i, ]
    indSocial <- SocialInfoMatrix[i, ]
    # create task vector to be output and bound
    taskThresh <- rep(0, length(indThresh))
    # loop through each task within individual
    for (j in 1:length(taskThresh)) {
      stim <- stimulusThisStep[j]
      thresh <- indThresh[j] - (epsilon * indSocial[j])
      if (stim > thresh) {
        taskThresh[j] <- 1
      }
    }
    return(taskThresh)
  })
  # bind and return
  thresholdP <- do.call("rbind", thresholdP)
  colnames(thresholdP) <- paste0("ThreshProb", 1:ncol(thresholdP))
  rownames(thresholdP) <- paste0("v-", 1:nrow(thresholdP))
  return(thresholdP)
}

####################
# Deterministic Threshold with social information
####################
adjustThresholdsSocial <- function(SocialNetwork, ThresholdMatrix, X_sub_g, epsilon) {
  # Calculate "sum" of task states/probs of neighbors
  NeighborSums <- t(SocialNetwork) %*% X_sub_g
  totalSums <- rowSums(NeighborSums)
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





