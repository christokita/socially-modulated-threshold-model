################################################################################
#
# Demonstrate sampling bias (for showing scaling of interaction bias with group size)
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")

####################
# Set variables
####################
beta <- 1.1 # interaction bias term
Ns <- c(10, 50, 100) # group sizes to test
gens <- 50000 #equivalnet of simulation run length

####################
# Set variables
####################
sample_results <- lapply(Ns, function(n) {
  # Establish "individuals:
  individuals <- data.frame(Individual = seq(1, n), 
                            Type = c(rep(1, n/2), rep(2, n/2)),
                            Weights = c(rep(beta, n/2), rep(1, n/2)))
  # Sample
  samples <- list()
  for (t in 1:gens) {
    if (t %% 10000 == 0) {
      print(t)
    }
    sampled_ind <- sample(individuals$Individual, 1, prob = individuals$Weights)
    row <- individuals[sampled_ind, 1:2]
    samples[[t]] <- row
  }
  samples <- do.call("rbind", samples)
  hist(samples$Individual)
})



