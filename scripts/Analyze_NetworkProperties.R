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
# Clustering coefficient
####################
g_example <- graph_from_adjacency_matrix(adjmatrix = soc_networks[[8]][[1]], weighted = T)
transitivity(graph = g_example, weights = E(g_example)$weight)



