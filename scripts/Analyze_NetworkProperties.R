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
# Set group size and replicate
size <- 35
size <- size/5
replicate <- 1

# Get graph
example_graph <- soc_networks[[size]][[replicate]]
example_thresh <- as.data.frame(thresh_data[[size]][[replicate]])
example_thresh$Id <- row.names(example_thresh)

# Get task distribution
example_task_dist <- task_dist[task_dist$n == size*5 & task_dist$replicate == replicate, ]
example_task_dist$total_activity <- rowSums(example_task_dist[ , 1:2])
mean_activity <- median(example_task_dist$total_activity)
example_task_dist$activity_level <- ifelse(example_task_dist$total_activity > mean_activity, "Above", "Below")

# Remove extra edges
percentiles <- quantile(example_graph, na.rm = TRUE)
fiftypercent <- percentiles[3]
seventyfivepercent <- percentiles[4]
example_graph[example_graph < fiftypercent] <- 0
diag(example_graph) <- 0

# Prep graph
g <- graph_from_adjacency_matrix(example_graph, mode = "undirected", weighted = TRUE)
edgelist <- get.edgelist(g)
edgelist <- as.data.frame(edgelist)
names(edgelist) <- c("Source", "Target")
edgelist$Weight <- E(g)$weight 
nodelist <- data.frame("Id" = row.names(example_graph))
nodelist$total_activity <- example_task_dist$total_activity
nodelist$activity_level <- example_task_dist$activity_level

# Output
write.csv(edgelist, file = paste0("output/Networks/ExampleNetworks/Activity_GroupSize", 5*size, ".csv"), row.names = FALSE)
write.csv(nodelist, file = paste0("output/Networks/ExampleNetworks/Activity_GroupSize", 5*size, "nodelist.csv"), row.names = FALSE)


# Plot with igraph
V(g)$color <- ifelse(example_task_dist$activity_level == "Above", "#d6604d", "#bababa")
graph_attr(g, "layout") <- layout_with_lgl
layout <- layout.forceatlas2(g)
plot(g, layout = layout)
