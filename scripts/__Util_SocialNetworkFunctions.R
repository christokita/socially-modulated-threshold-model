##################################################
#
# Social Network Functions
#
##################################################

####################
# Form connections probabilistically -  new formulation
####################
temporalNetwork <- function(X_sub_g, bias) {
  dimension <- nrow(X_sub_g)
  g_adj <- matrix(data = rep(0, dimension*dimension), ncol = dimension)
  # Loop through row individuals
  for (i in 1:dimension) {
    # get task performance of individual
    task <- which(X_sub_g[i, ] == 1)
    # set up list of potential connections
    potential <- seq(1:dimension)
    baseline_prob <- rep(1, length(potential))
    # loop through column individuals
    # If inactive, all connections equal prob
    if (length(task) == 0) {
      potential <- potential[-i] #remove self
      baseline_prob <- baseline_prob[-i] #remove self
      connection <- sample(x = potential, size = 1, prob = baseline_prob)
      g_adj[i, connection] <- 1
      g_adj[connection, i] <- 1
    } else { # If active, biased towards individuals in same state
      # find which individuals are perfoming same task and relatively weight probabilities
      same <- which(X_sub_g[ , task] == 1)
      baseline_prob[same] <- baseline_prob[same] * bias
      potential <- potential[-i] #remove self
      baseline_prob <- baseline_prob[-i] #remove self
      connection <- sample(x = potential, size = 1, prob = baseline_prob)
      g_adj[i, connection] <- 1
      g_adj[connection, i] <- 1
    }
  }
  # bind and name columns
  diag(g_adj) <- 0 #remove self-connections
  colnames(g_adj) <- paste0("v-", 1:dimension)
  rownames(g_adj) <- paste0("v-", 1:dimension)
  return(g_adj)
}


####################
# Form connections probabilistically
####################
temporalNetwork_OLD <- function(X_sub_g, p, bias) {
  dimension <- nrow(X_sub_g)
  g_adj <- matrix(data = rep(NA, dimension*dimension), ncol = dimension)
  # Loop through row individuals
  for (i in 1:dimension) {
    # get task performance of individual
    task <- which(X_sub_g[i, ] == 1)
    # loop through column individuals
    # If inactive, all connections equal prob
    if (length(task) == 0) {
      for (j in 1:dimension) {
        con <- sample(x = c(0, 1), size = 1, prob = c(1 - p, p))
        g_adj[i, j] <- con
        g_adj[j, i] <- con
      }
    }
    # If active, biased towards individuals in same state
    else {
      for (j in i:dimension) {
        same <- which(X_sub_g[ , task] == 1)
        # if same, biased
        if (j %in% same) {
          con <- sample(x = c(0, 1), size = 1, prob = c(1 - bias * p, bias * p))
          g_adj[i, j] <- con
          g_adj[j, i] <- con
        } else {
          con <- sample(x = c(0, 1), size = 1, prob = c(1 - p, p))
          g_adj[i, j] <- con
          g_adj[j, i] <- con
        }
      }
    }
  }
  # bind and name columns
  diag(g_adj) <- 0 #remove self-connections
  colnames(g_adj) <- paste0("v-", 1:dimension)
  rownames(g_adj) <- paste0("v-", 1:dimension)
  return(g_adj)
}

