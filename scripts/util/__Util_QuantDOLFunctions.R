##################################################
#
# Qunatifying Division of Labor Functions
#
##################################################

####################
# Mutual Entropy DOL Measure
####################
# From Gorelick, Bertram, Killeen, & Fewell (2004)
mutualEntropy <- function(TotalStateMat) {
  # Normalize matrix
  # normMat <- TotalStateMat / rowSums(TotalStateMat)
  normMat <- TotalStateMat / sum(TotalStateMat)
  # Total Individuals
  n <- nrow(normMat)
  m <- ncol(normMat)
  total <- sum(normMat)
  # Shannon's entropy of individuals H(X)
  H_x <- apply(normMat, MARGIN = 1, function(ind) {
    p_x <- sum(ind) / total
    h_x <- p_x * log2(p_x)
    if(is.na(h_x)) {
      h_x <- 0
    }
    return(h_x)
  })
  # Shannon's entropy of tasks H(Y)
  H_y <- apply(normMat, MARGIN = 2, function(task) {
    p_y <- sum(task) / total
    h_y <- p_y * log2(p_y)
    if(is.na(h_y)) {
      h_y <- 0
    }
    return(h_y)
  })
  # Mutual entropy I(X,Y)
  I_xy <- lapply(1:n, function(ind) {
    # Loop through tasks for each individual
    mutualEntr <- rep(NA, m)
    for (task in 1:m) {
      # joint probability p(x,y)
      p_xy <- normMat[ind, task] / total
      # calculate log portion
      p_x <- sum(normMat[ind, ]) / total
      p_y <- sum(normMat[ , task]) / total
      logVal <- log2(p_xy / (p_x * p_y))
      # If entry has zero probability, set total value to zero (instead of NA/-Inf)
      entry <- p_xy * logVal
      if (is.na(entry)) {
        entry <- 0
      }
      # enter into list
      mutualEntr[task] <- entry
    }
    mutualEntr <- sum(mutualEntr)
    return(mutualEntr)
  })
  # Sum values 
  H_x <- -sum(H_x)
  H_y <- -sum(H_y)
  I_xy <- sum(unlist(I_xy))
  # Calcualte symmetrid division of labor D(x,y)
  D_sym  <- I_xy / sqrt(H_x * H_y)
  D_task <- I_xy / H_x
  D_ind  <- I_xy / H_y
  # Dataframe
  D <- data.frame(Dsym = D_sym, Dtask = D_task, Dind = D_ind)
  D <- as.matrix(D)
  # Return 
  return(D)
}

####################
# Calcualte task rank
####################
calculateTaskRank <- function(TaskStepMat) {
  # Loop through columns
  for (column in 1:ncol(TaskStepMat)) {
    TaskStepMat[ , column] <- dense_rank(TaskStepMat[ , column])
  }
  # Return
  return(TaskStepMat)
}


