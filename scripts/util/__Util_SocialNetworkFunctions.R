##################################################
#
# Social Network Functions
#
##################################################

####################
# Form connections probabilistically -  new formulation
####################
temporalNetwork <- function(X_sub_g, prob_interact, bias) {
  dimension <- nrow(X_sub_g)
  g_adj <- matrix(data = rep(0, dimension*dimension), ncol = dimension)
  # Loop through row individuals
  for (i in 1:dimension) {
    # Determine if interact
    interact <- sample(x = c(1, 0), size = 1, prob = c(prob_interact, 1 - prob_interact))
    if (interact == 1) {
      # get task performance of individual
      task <- which(X_sub_g[i, ] == 1)
      # set up list of potential connections
      potential <- seq(1:dimension)
      baseline_prob <- rep(1, length(potential))
      # catch for if there is only two individuals
      if (length(baseline_prob) == 2) { 
        potential <- c(potential, 0)
        baseline_prob <- c(baseline_prob, 0)
      }
      # loop through column individuals
      # If inactive, all connections equal prob
      if (length(task) == 0) {
        potential <- potential[-i] # remove self
        baseline_prob <- baseline_prob[-i] # remove self
        connection <- sample(x = potential, size = 1, prob = baseline_prob)
        g_adj[i, connection] <- 1
        g_adj[connection, i] <- 1
      } else { # If active, biased towards individuals in same state
        # find which individuals are perfoming same task and relatively weight probabilities
        same <- which(X_sub_g[ , task] == 1)
        baseline_prob[same] <- baseline_prob[same] * bias
        potential <- potential[-i] # remove self
        baseline_prob <- baseline_prob[-i] # remove self
        connection <- sample(x = potential, size = 1, prob = baseline_prob)
        g_adj[i, connection] <- 1
        g_adj[connection, i] <- 1
      }
    }
  }
  # bind and name columns
  diag(g_adj) <- 0 #remove self-connections
  colnames(g_adj) <- paste0("v-", 1:dimension)
  rownames(g_adj) <- paste0("v-", 1:dimension)
  return(g_adj)
}


####################
# Form connections probabilistically -  new formulation
####################
temporalNetwork_OLD <- function(X_sub_g, bias) {
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
    #
    # FIX: fix baseline_prob vetor. technically fine, but reweight probs before feeding into sample. 
    #
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