################################################################################
#
# Analyze Networks of Social Interactions
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")
library(RColorBrewer)
library(scales)
library(igraph)

filename <- "Sigma0.0-Epsilon0.01-Bias1.1"

# Cutoff for threshold ratio to allow easier plotting
ThreshCutoffValue <- 10
ThreshCutoffReplacement <- Inf
ThreshCutoffReplacementColor <- 10

####################
# Load data
####################
# Load social
load("output/Rdata/Sigma0.0-Epsilon0.01-Bias1.1.Rdata")

soc_groups_graphs <- groups_graphs
soc_groups_threshMat <- groups_thresh
soc_graphs <- unlist(groups_graphs, recursive = FALSE)
soc_threshMat <- unlist(groups_thresh, recursive = FALSE)
soc_actMat <- unlist(groups_taskDist, recursive = FALSE)

social_graphs <- lapply(1:length(soc_graphs), function(i) {
  # Calculated degree
  degree <- rowSums(soc_graphs[[i]])
  degree <- as.data.frame(degree)
  degree$Id <- row.names(degree)
  # Calculate thresholds
  thresh <- soc_threshMat[[i]] %>% 
    as.data.frame(.) %>% 
    mutate(ThreshRatio = log(Thresh1 / Thresh2),
           Id = row.names(.)) %>% 
    mutate(ThreshRatioRaw = ThreshRatio)
  thresh$ThreshRatio[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffReplacement
  thresh$ThreshRatio[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffReplacement
  thresh$ThreshRatioColor <- thresh$ThreshRatio
  thresh$ThreshRatioColor[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffReplacementColor
  thresh$ThreshRatioColor[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffReplacementColor
  # Calculate actibity
  activity <- soc_actMat[[i]] %>% 
    as.data.frame(.) %>% 
    mutate(ActRatio = log(Task1 / Task2),
           ActTotal = Task1 + Task2,
           Id = row.names(.))
  # Merge and return
  mergedNodes <- merge(degree, thresh)
  mergedNodes <- merge(mergedNodes, activity)
  return(mergedNodes)
})

social_graphs <- do.call("rbind", social_graphs)

# Load fixed
load("output/Rdata/Sigma0.05-FIXED-Bias1.1.Rdata")

fix_groups_graphs <- groups_graphs
fix_groups_threshMat <- groups_thresh
fix_graphs <- unlist(groups_graphs, recursive = FALSE)
fix_threshMat <- unlist(groups_thresh, recursive = FALSE)
fix_actMat <- unlist(groups_taskDist, recursive = FALSE)

fixed_graphs <- lapply(1:length(fix_graphs), function(i) {
  # Calculated degree
  degree <- rowSums(fix_graphs[[i]])
  degree <- as.data.frame(degree)
  degree$Id <- row.names(degree)
  # Calculate thresholds
  thresh <- fix_threshMat[[i]] %>% 
    as.data.frame(.) %>% 
    mutate(ThreshRatio = log(Thresh1 / Thresh2),
           Id = row.names(.)) %>% 
    mutate(ThreshRatioRaw = ThreshRatio)
  thresh$ThreshRatio[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffReplacement
  thresh$ThreshRatio[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffReplacement
  thresh$ThreshRatioColor <- thresh$ThreshRatio
  thresh$ThreshRatioColor[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffReplacementColor
  thresh$ThreshRatioColor[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffReplacementColor
  # Calculate actibity
  activity <- fix_actMat[[i]] %>% 
    as.data.frame(.) %>% 
    mutate(ActRatio = log(Task1 / Task2),
           ActTotal = Task1 + Task2,
           Id = row.names(.))
  # Merge and return
  mergedNodes <- merge(degree, thresh)
  mergedNodes <- merge(mergedNodes, activity)
  return(mergedNodes)
})

fixed_graphs <- do.call("rbind", fixed_graphs)


###### Interaction Matrix Plot ###### 

type_groups_graphs <- fix_groups_graphs
type_groups_threshMat <- fix_groups_threshMat

interaction_graphs <- lapply(1:length(type_groups_graphs), function(i) {
  # Get graphs
  graphs <- type_groups_graphs[[i]]
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
    # this_graph <- scale(this_graph) # normalize relative to the average
    not_chosen <- 1 - ( 1 / (dimensions[1] - 1) )
    expected_random <-  1 - not_chosen^2
    this_graph <- (this_graph - expected_random) / expected_random #relative to expected by random (i.e., 1 - chance of not being chosen^2)
    this_graph <- matrix(data = this_graph, nrow = dimensions[1], ncol = dimensions[2])
    colnames(this_graph) <- labs
    rownames(this_graph) <- labs
    # Calculate thresh ratio
    # ind <- replicates * i + j
    thresh <- as.data.frame(type_groups_threshMat[[i]][j])
    thresh$ThreshRatio <- log(thresh$Thresh1 / thresh$Thresh2)
    ratio <- order(thresh$ThreshRatio)
    # Create order by threshold ratio
    this_graph <- this_graph[ratio, ratio]
    colnames(this_graph) <- paste0("t-", 1:dimensions[1])
    rownames(this_graph) <- paste0("t-", 1:dimensions[1])
    # return
    return(this_graph)
  })
  # Avearge across all to make 'typical' adjacency matrix
  avg_g <- Reduce("+", size_graph) / length(size_graph)
  # Create graph object
  g <- graph_from_adjacency_matrix(avg_g, mode = c("directed"), weighted = TRUE, diag = TRUE)
  # Get node and edge list
  node_list <- get.data.frame(g, what = "vertices")
  edge_list <- get.data.frame(g, what = "edges")
  # Create dataframe for plotting
  plot_data <- edge_list %>% mutate(
    to = factor(to, levels = node_list$name),
    from = factor(from, levels = node_list$name))
  # Plot
  groupsize <- ncol(avg_g)
  gg_avg_adj <- ggplot(plot_data, aes(x = rev(from), y = to, fill = weight)) +
    geom_raster() +
    theme_bw() +
    # Because we need the x and y axis to display every node,
    # not just the nodes that have connections to each other,
    # make sure that ggplot does not drop unused factor levels
    scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
    # scale_fill_gradientn(colours = rev(brewer.pal(9,"RdYlBu")), na.value = "white", limit = c(-1.5, 1.5), oob = squish) +
    scale_fill_gradientn(name = "Relative Interaction\nFrequency",
                         # colours = c(brewer.pal(7,"PiYG")),
                         colours = c('#525252','#5b5b5b','#646464','#6e6e6e','#787878','#818181','#8b8b8b','#959595','#a0a0a0','#a9a9a9','#b4b4b4','#bfbfbf','#c8c8c8','#d4d4d4','#dedede','#e9e9e9','#f4f4f4','#ffffff','#edf5f9','#dee9f2','#d3ddec','#c7d1e5','#bfc4de','#b7b7d7','#b0aad0','#a99ec8','#a391c1','#9e83b9','#9a76b1','#9569a9','#915aa1','#8c4c98','#893c8f','#852986','#810f7c'),
                         # colours = rev(c("#F6BDAA", "#EC8591", "#E15287", "#AC3987", "#6B249C", "#4D1B7A", "#381B4A")),
                         na.value = "black", 
                         limit = c(-0.05, 0.05),
                         # limit = c(0.95, 1.05),
                         oob = squish) +
    theme(axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      aspect.ratio = 1,
      # Hide the legend (optional)
      legend.position = "none",
      legend.key.height = unit(0.38, "in"),
      panel.border = element_rect(size = 1.5),
      title = element_blank()) +
    ggtitle(paste0("Group Size = ", groupsize))
  # return graph
  return(gg_avg_adj)
  # return(avg_g)
})

png(filename = paste0("output/Networks/IntMat_", filename, ".png"), width = 5, height = 2.5, units = "in", res = 600)
multiplot(plotlist = interaction_graphs, layout = matrix(c(seq(1:length(interaction_graphs))), nrow=2, byrow=TRUE))
dev.off()

# Subset of graphs plot for figure
png(filename = paste0("output/Networks/FIGURE_IntMat_Single_", filename, ".png"), width = 2, height = 2, units = "in", res = 600)
multiplot(plotlist = interaction_graphs[c(7)], layout = matrix(c(seq(1:length(interaction_graphs[c(1, 2, 5, 8)]))), nrow = 1, byrow=TRUE))
dev.off()

# Subset of graphs plot for figure
png(filename = paste0("output/Networks/FIGURE_IntMat_", filename, ".png"), width = 4, height = 1, units = "in", res = 600)
multiplot(plotlist = interaction_graphs[c(2, 4, 6, 8)], layout = matrix(c(seq(1:length(interaction_graphs[c(1, 2, 5, 8)]))), nrow = 1, byrow=TRUE))
dev.off()

