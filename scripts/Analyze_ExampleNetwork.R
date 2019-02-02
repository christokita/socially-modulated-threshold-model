################################################################################
#
# Analyzing single example network
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

p <- 1 #prob of interact
run <- "Sigma0-Epsilon0.1-Beta1.1"

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

####################
# Output example graph
####################
# Set threshold max/min
thresh_limit <- 100
# Set group size and replicate
size <- 80
size <- size/5
replicate <- 1
# Get graph
example_graph <- soc_networks[[size]][[replicate]]
example_thresh <- as.data.frame(thresh_data[[size]][[replicate]])
example_thresh$Id <- row.names(example_thresh)
example_thresh$ThreshBias <- example_thresh$Thresh1 - example_thresh$Thresh2
example_thresh$ThreshBiasBounded <- example_thresh$ThreshBias
example_thresh$ThreshBiasBounded[example_thresh$ThreshBiasBounded < -thresh_limit] <- -thresh_limit
example_thresh$ThreshBiasBounded[example_thresh$ThreshBiasBounded > thresh_limit] <- thresh_limit
# If no node reaches upper or lower limits, add for coloring purposes in gephi
if (sum(example_thresh$ThreshBias == thresh_limit) == 0) {
  max_row <- data.frame(Thresh1 = 50, Thresh2 = 50,
                        n = size * 5, sim = 0, chunk = 0,
                        Id = "Max", ThreshBias = thresh_limit, ThreshBiasBounded = thresh_limit)
  example_thresh <- rbind(example_thresh, max_row)
}
if (sum(example_thresh$ThreshBias == -thresh_limit) == 0) {
  min_row <- data.frame(Thresh1 = 50, Thresh2 = 50,
                        n = size * 5, sim = 0, chunk = 0,
                        Id = "Min", ThreshBias = -thresh_limit, ThreshBiasBounded = -thresh_limit)
  example_thresh <- rbind(example_thresh, min_row)
}
# example_thresh$ThreshRatio <- log(example_thresh$Thresh1 / example_thresh$Thresh2)
# example_thresh$ThreshRatioBounded <- example_thresh$ThreshRatio
# example_thresh$ThreshRatioBounded[example_thresh$ThreshRatioBounded < -thresh_limit] <- -thresh_limit
# example_thresh$ThreshRatioBounded[example_thresh$ThreshRatioBounded > thresh_limit] <- thresh_limit
# # If no node reaches upper or lower limits, add for coloring purposes in gephi
# if (sum(example_thresh$ThreshRatioBounded == thresh_limit) == 0) {
#   max_row <- data.frame(Thresh1 = NA, Thresh2 = NA, 
#                         n = NA, sim = NA, chunk = NA,
#                         Id = "Max", ThreshRatio = thresh_limit, ThreshRatioBounded = thresh_limit)
#   example_thresh <- rbind(example_thresh, max_row)
# }
# if (sum(example_thresh$ThreshRatioBounded == -thresh_limit) == 0) {
#   min_row <- data.frame(Thresh1 = NA, Thresh2 = NA, 
#                         n = NA, sim = NA, chunk = NA,
#                         Id = "Min", ThreshRatio = -thresh_limit, ThreshRatioBounded = -thresh_limit)
#   example_thresh <- rbind(example_thresh, min_row)
# }
# # Calculate values expected
# Zero out interactions equal to or less than random
not_chosen <- 1 - (( 1 / (nrow(example_graph) - 1)) * p)
expected_random <-  1 - not_chosen^2
example_graph[example_graph <= expected_random] <- 0

# Or take in those in top X percentile
percentiles <- quantile(example_graph, na.rm = TRUE)
fiftypercent <- percentiles[3]
seventyfivepercent <- percentiles[4]
example_graph[example_graph <= fiftypercent] <- 0
diag(example_graph) <- 0

# Turn into graph object to get edgelist
g <- graph_from_adjacency_matrix(example_graph, mode = "undirected", weighted = TRUE)
edgelist <- get.edgelist(g)
edgelist <- as.data.frame(edgelist)
names(edgelist) <- c("Source", "Target")
edgelist$Weight <- E(g)$weight 

# Write
write.csv(edgelist, file = paste0("output/Networks/ExampleNetworks/GroupSize", 5*size, "_AboveRandom_edgelist.csv"), row.names = FALSE)
write.csv(example_thresh, file = paste0("output/Networks/ExampleNetworks/GroupSize", 5*size, "nodelist.csv"), row.names = FALSE)

