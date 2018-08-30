################################################################################
#
# Analyzing social network features comapred to activity and other factors
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

p <- 1 #prob of interact
run <- "Sigma0-Epsilon0.1-Beta1.1"

####################
# Load and process data
####################
# Load social networks
files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", run, "/"), full.names = TRUE)
soc_networks <- list()
for (i in 1:length(files)) {
  load(files[i])
  soc_networks[[i]] <- listed_data
}

# Load threshold matrices
files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh/", run, "/"), full.names = TRUE)
thresh_data <- list()
for (i in 1:length(files)) {
  load(files[i])
  thresh_data[[i]] <- listed_data
}

# Load activity profiles
load(paste0("output/Rdata/_ProcessedData/TaskDist/", run, ".Rdata"))
task_dist <- compiled_data
task_dist$replicate <- task_dist$sim * task_dist$chunk
rm(compiled_data)


####################
# Total interactions vs. activity
####################
# Create space for data collection
task_activ <- task_dist
task_activ$degree <- NA

# Get data from social networks
for(i in 1:length(soc_networks)) {
  # Get group size graphs
  graphs <- soc_networks[[i]]
  # Get individual graphs
  for(j in 1:length(graphs)) {
    this_graph <- graphs[[j]]
    degrees <- rowSums(this_graph)
    # add to dataframe
    n <- 5 * i
    task_activ$degree[task_activ$n == n & task_activ$replicate == j] <- degrees
  }
}

# Plot
gg_act_net <- ggplot(data = task_activ, aes(x = Task1, y = degree, colour = n)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~n)
gg_act_net


####################
# Network plot vs. above/below average activity
####################
interaction_graphs <- lapply(1:length(soc_networks), function(i) {
  # Get graphs
  graphs <- soc_networks[[i]]
  replicates <- length(graphs)
  # For each each compute interaction matrix
  # Get graph and make adjacency matrix
  size_graph <- lapply(1:length(graphs), function(j) {
    # Format: set diagonal, rescale, and make adj matrix
    this_graph <- graphs[[j]]
    diag(this_graph) <- NA
    dimensions <- dim(this_graph)
    labs <- colnames(this_graph)
    this_graph <- as.vector(this_graph)
    not_chosen <- 1 - (( 1 / (dimensions[1] - 1)) * p)
    expected_random <-  1 - not_chosen^2
    this_graph <- (this_graph - expected_random) / expected_random #relative to expected by random (i.e., 1 - chance of not being chosen^2)
    this_graph <- matrix(data = this_graph, nrow = dimensions[1], ncol = dimensions[2])
    colnames(this_graph) <- labs
    rownames(this_graph) <- labs
    # Calculate thresh ratio
    # ind <- replicates * i + j
    thresh <- as.data.frame(thresh_data[[i]][j])
    thresh$ThreshRatio <- log(thresh$Thresh1 / thresh$Thresh2)
    ratio <- order(thresh$ThreshRatio)
  })
})

# Plot
gg_act_net <- ggplot(data = task_activ, aes(x = Task1, y = degree, colour = n)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~n)
gg_act_net

