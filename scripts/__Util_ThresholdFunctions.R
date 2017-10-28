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
# Threshold function
####################
threshProb <- function(s, phi, nSlope) {
  T_vi <- (s^nSlope) / (s^nSlope + phi^nSlope)
}


####################
# Output Threshold Demands
####################
calcThresholdProbMat <- function(TimeStep, ThresholdMatrix, StimulusMatrix, nSlope) {
  # select proper stimulus for this time step
  stimulusThisStep <- StimulusMatrix[TimeStep, ]
  # calculate threshold probabilities for one individual
  thresholdP <- lapply(1:nrow(ThresholdMatrix), function(i) {
    # select row for individual in threshold matrix
    indThresh <- ThresholdMatrix[i, ]
    # create task vector to be output and bound
    taskThresh <- rep(NA, length(indThresh))
    # loop through each task within individual
    for (j in 1:length(taskThresh)) {
      taskThresh[j] <- threshProb(s = stimulusThisStep[j], phi = indThresh[j], nSlope = nSlope)
    }
    return(taskThresh)
  })
  # bind and return
  thresholdP <- do.call("rbind", thresholdP)
  thresholdP[is.na(thresholdP)] <- 0 #fix NAs where thresh was 0 and stim was 0
  colnames(thresholdP) <- paste0("ThreshProb", 1:ncol(thresholdP))
  rownames(thresholdP) <- paste0("v-", 1:nrow(thresholdP))
  return(thresholdP)
}

####################
# Output Threshold Demands
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
# Exhaustion Threshold Function
####################
calcExhaustDemand <- function(ExhaustStim, ExhaustThreshVector, nSlope) {
  # Create output matrix
  exhaustUpdate <- matrix(rep(NA, length(ExhaustStim)), ncol = 1)
  colnames(exhaustUpdate) <- "ExhaustProb"
  rownames(exhaustUpdate) <- paste0("v-", 1:length(ExhaustStim))
  # Loop through individuals
  for (i in 1:length(exhaustUpdate)) {
    exhaustUpdate[i] <- threshProb(s = ExhaustStim[i], phi = ExhaustThreshVector[i], nSlope = nSlope)
  }
  # Return
  return(exhaustUpdate)
}

####################
# Self-reinforcement of Threshold
####################
adjustThresholds <- function(ThresholdMatrix, X_sub_g, phi, gamma, lowerThresh, upperThresh) {
  for (i in 1:nrow(X_sub_g)) {
    for (j in 1:ncol(X_sub_g)) {
      adjust <- ((1 - X_sub_g[i, j]) * phi - X_sub_g[i, j] * gamma)
      ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
      if (ThresholdMatrix[i, j] < lowerThresh) {
        ThresholdMatrix[i, j] <- lowerThresh
      } else if (ThresholdMatrix[i, j] > upperThresh) {
        ThresholdMatrix[i, j] <- upperThresh
      }
    }
  }
  return(ThresholdMatrix)
}

####################
# Self-reinforcement of Threshold
####################
adjustThresholdsSocial <- function(ThresholdMatrix, X_sub_g, phi, gamma, c, d, SocialNetwork, lowerThresh, upperThresh, Average) {
  # Calculate "sum" of task states/probs of neighbors
  NeighborSums <- SocialNetwork %*% X_sub_g
  # Calculate total neighbors (degree of individual)
  DegSum <- rowSums(SocialNetwork)
  # Calculate frequency performing
  freqPerforming <- NeighborSums / DegSum
  freqPerforming[is.na(freqPerforming)] <- 0
  # Loop through individuals
  # If summing interactions
  if (Average == FALSE) {
    for (i in 1:nrow(X_sub_g)) {
      for (j in 1:ncol(X_sub_g)) {
        # Establish variables
        state <-  X_sub_g[i, j]
        sum <- NeighborSums[i, j]
        # Adjust
        adjust <- ( ((1 - state) * phi * (1 - d * sum)) - (state * gamma * (1 + c * sum)) )
        ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
        if (ThresholdMatrix[i, j] < lowerThresh) {
          ThresholdMatrix[i, j] <- lowerThresh
        } else if (ThresholdMatrix[i, j] > upperThresh) {
          ThresholdMatrix[i, j] <- upperThresh
        }
      }
    }
  } else if (Average == TRUE) {
    for (i in 1:nrow(X_sub_g)) {
      for (j in 1:ncol(X_sub_g)) {
        # Establish variables
        state <-  X_sub_g[i, j]
        freq <- freqPerforming[i, j]
        # Adjust
        adjust <- ( ((1 - state) * phi * (1 - d * freq)) - (state * gamma * (1 + c * freq)) )
        ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
        if (ThresholdMatrix[i, j] < lowerThresh) {
          ThresholdMatrix[i, j] <- lowerThresh
        } else if (ThresholdMatrix[i, j] > upperThresh) {
          ThresholdMatrix[i, j] <- upperThresh
        }
      }
    }
  }
  return(ThresholdMatrix)
}

####################
# Self-reinforcement of Threshold
####################
adjustThresholdsSocial_bound <- function(ThresholdMatrix, X_sub_g, phi, gamma, c, d, SocialNetwork, lowerThresh, upperThresh, Average) {
  # Calculate "sum" of task states/probs of neighbors
  NeighborSums <- SocialNetwork %*% X_sub_g
  # Calculate total neighbors (degree of individual)
  DegSum <- rowSums(SocialNetwork)
  # Calculate frequency performing
  freqPerforming <- NeighborSums / DegSum
  freqPerforming[is.na(freqPerforming)] <- 0
  # Loop through individuals
  # If summing interactions
  if (Average == FALSE) {
    for (i in 1:nrow(X_sub_g)) {
      for (j in 1:ncol(X_sub_g)) {
        # Establish variables
        state <-  X_sub_g[i, j]
        sum <- NeighborSums[i, j]
        # Adjust
        adjust <- ( ((1 - state) * phi * (1 - d * sum)) - (state * gamma * (1 + c * sum)) )
        ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
        if (ThresholdMatrix[i, j] < lowerThresh[i, j]) {
          ThresholdMatrix[i, j] <- lowerThresh[i, j]
        } else if (ThresholdMatrix[i, j] > upperThresh[i, j]) {
          ThresholdMatrix[i, j] <- upperThresh[i, j]
        }
      }
    }
  } else if (Average == TRUE) {
    for (i in 1:nrow(X_sub_g)) {
      for (j in 1:ncol(X_sub_g)) {
        # Establish variables
        state <-  X_sub_g[i, j]
        freq <- freqPerforming[i, j]
        # Adjust
        adjust <- ( ((1 - state) * phi * (1 - d * freq)) - (state * gamma * (1 + c * freq)) )
        ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
        if (ThresholdMatrix[i, j] < lowerThresh[i, j]) {
          ThresholdMatrix[i, j] <- lowerThresh[i, j]
        } else if (ThresholdMatrix[i, j] > upperThresh[i, j]) {
          ThresholdMatrix[i, j] <- upperThresh[i, j]
        }
      }
    }
  }
  return(ThresholdMatrix)
}

####################
# Self-reinforcement of Threshold
####################
adjustThresholdsSocialInhibition <- function(ThresholdMatrix, X_sub_g, phi, gamma, c, d, SocialNetwork, lowerThresh, upperThresh, Average) {
  # Calculate "sum" of task states/probs of neighbors
  NeighborSums <- SocialNetwork %*% X_sub_g
  # Calculate total neighbors (degree of individual)
  DegSum <- rowSums(SocialNetwork)
  # Calculate Inactive
  Inactive <- DegSum - rowSums(NeighborSums)
  InactiveFreq <- Inactive / DegSum
  InactiveFreq[is.na(InactiveFreq)] <- 0
  # Calculate frequency performing
  freqPerforming <- NeighborSums / DegSum
  freqPerforming[is.na(freqPerforming)] <- 0
  # Loop through individuals
  # If summing interactions
  if (Average == FALSE) {
    for (i in 1:nrow(X_sub_g)) {
      for (j in 1:ncol(X_sub_g)) {
        # Establish variables
        state <-  X_sub_g[i, j]
        absDiff <- NeighborSums[i, j] - (DegSum[i] - Inactive[i] - NeighborSums[i, j])
        # Adjust
        adjust <- ( ((1 - state) * phi * (1 - d * absDiff)) - (state * gamma * (1 + c * absDiff)) )
        ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
        if (ThresholdMatrix[i, j] < lowerThresh) {
          ThresholdMatrix[i, j] <- lowerThresh
        } else if (ThresholdMatrix[i, j] > upperThresh) {
          ThresholdMatrix[i, j] <- upperThresh
        }
      }
    }
  } else if (Average == TRUE) {
    for (i in 1:nrow(X_sub_g)) {
      for (j in 1:ncol(X_sub_g)) {
        # Establish variables
        state <-  X_sub_g[i, j]
        freqDiff <- freqPerforming[i, j] - (1 - InactiveFreq[i] - freqPerforming[i, j])
        # Adjust
        adjust <- ( ((1 - state) * phi * (1 - d * freqDiff)) - (state * gamma * (1 + c * freqDiff)) )
        ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
        if (ThresholdMatrix[i, j] < lowerThresh) {
          ThresholdMatrix[i, j] <- lowerThresh
        } else if (ThresholdMatrix[i, j] > upperThresh) {
          ThresholdMatrix[i, j] <- upperThresh
        }
      }
    }
  }
  return(ThresholdMatrix)
}


####################
# Plastic bound self-reinforcement of Threshold
####################
adjustThresholdsSocialInhibition_Bound <- function(ThresholdMatrix, X_sub_g, phi, gamma, c, d, SocialNetwork, lowerThresh, upperThresh, Average) {
  # Calculate "sum" of task states/probs of neighbors
  NeighborSums <- SocialNetwork %*% X_sub_g
  # Calculate total neighbors (degree of individual)
  DegSum <- rowSums(SocialNetwork)
  # Calculate Inactive
  Inactive <- DegSum - rowSums(NeighborSums)
  InactiveFreq <- Inactive / DegSum
  InactiveFreq[is.na(InactiveFreq)] <- 0
  # Calculate frequency performing
  freqPerforming <- NeighborSums / DegSum
  freqPerforming[is.na(freqPerforming)] <- 0
  # Loop through individuals
  # If summing interactions
  if (Average == FALSE) {
    for (i in 1:nrow(X_sub_g)) {
      for (j in 1:ncol(X_sub_g)) {
        # Establish variables
        state <-  X_sub_g[i, j]
        absDiff <- NeighborSums[i, j] - (DegSum[i] - Inactive[i] - NeighborSums[i, j])
        # Adjust
        adjust <- ( ((1 - state) * phi * (1 - d * absDiff)) - (state * gamma * (1 + c * absDiff)) )
        ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
        if (ThresholdMatrix[i, j] < lowerThresh[i, j]) {
          ThresholdMatrix[i, j] <- lowerThresh[i, j]
        } else if (ThresholdMatrix[i, j] > upperThresh[i, j]) {
          ThresholdMatrix[i, j] <- upperThresh[i, j]
        }
      }
    }
  } else if (Average == TRUE) {
    for (i in 1:nrow(X_sub_g)) {
      for (j in 1:ncol(X_sub_g)) {
        # Establish variables
        state <-  X_sub_g[i, j]
        freqDiff <- freqPerforming[i, j] - (1 - InactiveFreq[i] - freqPerforming[i, j])
        # Adjust
        adjust <- ( ((1 - state) * phi * (1 - d * freqDiff)) - (state * gamma * (1 + c * freqDiff)) )
        ThresholdMatrix[i, j] <- ThresholdMatrix[i, j] + adjust
        if (ThresholdMatrix[i, j] < lowerThresh[i, j]) {
          ThresholdMatrix[i, j] <- lowerThresh[i, j]
        } else if (ThresholdMatrix[i, j] > upperThresh[i, j]) {
          ThresholdMatrix[i, j] <- upperThresh[i, j]
        }
      }
    }
  }
  return(ThresholdMatrix)
}




