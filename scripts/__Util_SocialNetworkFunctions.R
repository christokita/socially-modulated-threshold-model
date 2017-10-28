##################################################
#
# Social Network Functions
#
##################################################

####################
# Form connections probabilistically
####################
temporalNetwork <- function(X_sub_g, p, bias) {
  dimension <- nrow(X_sub_g)
  g_adj <- matrix(data = rep(NA, dimension*dimension), ncol = dimension)
  # Loop through row individuals
  for (i in 1:dimension) {
    # get task performance of individual
    task <- which(X_sub_g[i, ] == 1)
    # loop through column individuals
    # If inactive, all connections equal prob
    if (length(task) == 0) {
      for (j in i:dimension) {
        con <- sample(x = c(0, 1), size = 1, prob = c(1 - p, p))
        g_adj[i, j] <- con
        g_adj[j, i] <- con
      }
    }
    # If active, biased towards individuals in same state
    else {
      for (j in 1:dimension) {
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

