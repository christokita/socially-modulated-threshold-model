################################################################################
#
# Analyzing the social networks of simulations
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


####################
# Graph relative interactrion rates
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
    # ratio <- order(thresh$ThreshRatio)
    ratio <- order(thresh$Thresh1)
    # Create order by threshold ratio
    this_graph <- this_graph[ratio, ratio]
    colnames(this_graph) <- paste0("v-", 1:dimensions[1])
    rownames(this_graph) <- paste0("v-", 1:dimensions[1])
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
  # Get info for plot
  groupsize <- ncol(avg_g)
  # if (groupsize < 20) {
  #   breaks <- c(1, seq(5, length(plot_data$to), 5))
  # } else if(groupsize < 35) {
  #   breaks <- c(1, seq(10, length(plot_data$to), 15))
  # } else {
  #   breaks <- c(1, seq(20, length(plot_data$to), 20))
  # }
  breaks <- c(1, length(unique(plot_data$to)))
  # Plot
  gg_avg_adj <- ggplot(plot_data, aes(x = from, y = to, fill = weight, color = weight)) +
    geom_tile() +
    theme_bw() +
    # Because we need the x and y axis to display every node,
    # not just the nodes that have connections to each other,
    # make sure that ggplot does not drop unused factor levels
    scale_x_discrete(drop = FALSE, expand = c(0, 0), 
                     position = "top",
                     breaks = levels(plot_data$to)[breaks]) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0), 
                     limits = rev(levels(plot_data$to)),
                     breaks = levels(plot_data$to)[breaks]) +
    scale_fill_gradientn(name = "Relative\ninteraction\nfrequency",
                         colours = c('#525252','#5b5b5b','#646464','#6e6e6e','#787878','#818181','#8b8b8b','#959595','#a0a0a0','#a9a9a9','#b4b4b4','#bfbfbf','#c8c8c8','#d4d4d4','#dedede','#e9e9e9','#f4f4f4','#ffffff','#edf5f9','#dee9f2','#d3ddec','#c7d1e5','#bfc4de','#b7b7d7','#b0aad0','#a99ec8','#a391c1','#9e83b9','#9a76b1','#9569a9','#915aa1','#8c4c98','#893c8f','#852986','#810f7c'),
                         # colours = rev(c("#F6BDAA", "#EC8591", "#E15287", "#AC3987", "#6B249C", "#4D1B7A", "#381B4A")),
                         na.value = "white", 
                         limit = c(-0.05, 0.05),
                         # limit = c(0.95, 1.05),
                         oob = squish) +
    scale_color_gradientn(name = "Relative\ninteraction\nfrequency",
                         colours = c('#525252','#5b5b5b','#646464','#6e6e6e','#787878','#818181','#8b8b8b','#959595','#a0a0a0','#a9a9a9','#b4b4b4','#bfbfbf','#c8c8c8','#d4d4d4','#dedede','#e9e9e9','#f4f4f4','#ffffff','#edf5f9','#dee9f2','#d3ddec','#c7d1e5','#bfc4de','#b7b7d7','#b0aad0','#a99ec8','#a391c1','#9e83b9','#9a76b1','#9569a9','#915aa1','#8c4c98','#893c8f','#852986','#810f7c'),
                         # colours = rev(c("#F6BDAA", "#EC8591", "#E15287", "#AC3987", "#6B249C", "#4D1B7A", "#381B4A")),
                         na.value = "white", 
                         limit = c(-0.05, 0.05),
                         # limit = c(0.95, 1.05),
                         oob = squish) +
    xlab("Individual") +
    ylab("Individual") +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          # axis.text = element_text(colour = "black", size = 6),
          # axis.title = element_text(size = 7),
          # axis.ticks = element_line(size = 0.3),
          aspect.ratio = 1,
          # Hide the legend (optional)
          legend.position = "right",
          legend.key.width = unit(3, "mm"),
          legend.key.height = unit(6, "mm"),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          # panel.background = element_rect(size = 0.3, fill = NA),
          panel.border = element_rect(size = 0.3, fill = NA),
          plot.title = element_blank()) +
    ggtitle(paste0("Group Size = ", groupsize))
  # return graph
  return(gg_avg_adj)
  # return(avg_g)
})

# Save
plots <- seq(4, 7, 1)
for (plot in plots) {
  gg_inter <- interaction_graphs[[plot]]
  ggsave(gg_inter, 
         filename = paste0("output/Networks/RawPlots/", run, "_n", plot*5 ,".svg"), 
         width = 20, 
         height = 20,
         units = "mm")
}

# save with legend
gg_inter <- interaction_graphs[1]
gg_inter
ggsave("output/Networks/RawPlots/LegendPlot.svg", width = 45, height = 45, units = "mm")

####################
# Graph relative interactrion rates: simplified rate
####################
simple_graphs <- lapply(1:length(soc_networks), function(i) {
  # Get graphs
  graphs <- soc_networks[[i]]
  replicates <- length(graphs)
  # For each each compute interaction matrix
  # Get graph and make adjacency matrix
  size_graph <- lapply(1:length(graphs), function(j) {
    # Format: set diagonal, rescale, and make adj matrix
    this_graph <- graphs[[j]]
    diag(this_graph) <- NA
    thresh <- as.data.frame(thresh_data[[i]][j])
    thresh$ThreshRatio <- log(thresh$Thresh1 / thresh$Thresh2)
    ratio <- order(thresh$ThreshRatio)
    # Create order by threshold ratio
    this_graph <- this_graph[ratio, ratio]
    colnames(this_graph) <- 1:nrow(this_graph)
    rownames(this_graph) <- colnames(this_graph)
    g <- graph.adjacency(adjmatrix = this_graph, weighted = T)
    edgelist_graph <- as.data.frame(get.edgelist(g))
    names(edgelist_graph) <- c("From", "To")
    edgelist_graph$Weight <- E(g)$weight 
    edgelist_graph$Interaction <- paste0(edgelist_graph$From, "-", edgelist_graph$To)
    # return
    return(edgelist_graph)
  })
  #Calculate baseline probability of interaction
  dimensions <- dim(graphs[[1]])
  not_chosen <- 1 - (( 1 / (dimensions[1] - 1)) * p)
  expected_random <-  1 - not_chosen^2
  # Bind
  all_edgelist <- do.call("rbind", size_graph)
  # Get pvalues of interaction rate relative to expected by random
  nonNA_edgelist <- all_edgelist %>% 
    group_by(Interaction) %>% 
    filter(!is.na(Weight)) %>% 
    mutate(pvalue = t.test(x = Weight, mu = expected_random)[3]) %>% 
    mutate(Significant = pvalue < 0.05) %>% 
    filter(Significant == TRUE)
  # Get the CI for those that are significantly different and then determine if greater or less than expected
  if (nrow(nonNA_edgelist) > 0) {
    sigDiff_edgelist <- nonNA_edgelist %>% 
    group_by(Interaction) %>% 
    mutate(CImin = t.test(x = Weight, mu = expected_random)[[4]][1] - expected_random,
           CImax = t.test(x = Weight, mu = expected_random)[[4]][1] - expected_random) %>% 
    mutate(DiffDirection = ifelse(test = CImin > 0 & CImax > 0, yes = 1, no = -1)) %>% 
    select(Interaction, DiffDirection) %>% 
    unique(.) 
    calc_edgelist <- merge(all_edgelist, sigDiff_edgelist, all = TRUE)
  } else {
    calc_edgelist <- all_edgelist
    calc_edgelist$DiffDirection <- NA
  }
  # Merge back together and then form graph
  calc_edgelist$DiffDirection[is.na(calc_edgelist$DiffDirection) & !is.na(calc_edgelist$Weight)] <- 0
  calc_edgelist <- calc_edgelist %>% 
    select(From, To, Interaction, DiffDirection) %>% 
    unique(.) %>% 
    mutate(From = as.numeric(as.character(From)),
           To = as.numeric(as.character(To))) %>% 
    arrange(From, To)
  avg_g <- matrix(data = calc_edgelist$DiffDirection, nrow = dimensions[1], byrow = T)
  colnames(avg_g) <- paste0("t-", 1:dimensions[1])
  rownames(avg_g) <- paste0("t-", 1:dimensions[1])
  # Create graph object
  g <- graph_from_adjacency_matrix(avg_g, mode = c("directed"), weighted = TRUE, diag = TRUE)
  # Get node and edge list
  node_list <- get.data.frame(g, what = "vertices")
  edge_list <- get.data.frame(g, what = "edges")
  # Create dataframe for plotting
  plot_data <- edge_list %>% mutate(
    to = factor(to, levels = rev(node_list$name)),
    from = factor(from, levels = node_list$name))
  # Get info for plot
  groupsize <- ncol(avg_g)
  if (groupsize < 30) {
    breaks <- c(1, seq(5, length(plot_data$to), 5))
  } else if(groupsize < 55) {
    breaks <- c(1, seq(10, length(plot_data$to), 10))
  } else {
    breaks <- c(1, seq(20, length(plot_data$to), 20))
  }
  # Plot
  gg_avg_adj <- ggplot(plot_data, aes(x = from, y = to, fill = weight)) +
    geom_tile() +
    theme_bw() +
    # Because we need the x and y axis to display every node,
    # not just the nodes that have connections to each other,
    # make sure that ggplot does not drop unused factor levels
    scale_x_discrete(drop = FALSE, expand = c(0, 0), 
                     position = "top",
                     breaks = levels(plot_data$to)[breaks]) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0), 
                     limits = levels(plot_data$to),
                     breaks = levels(plot_data$to)[breaks]) +
    # scale_fill_gradientn(colours = rev(brewer.pal(9,"RdYlBu")), na.value = "white", limit = c(-1.5, 1.5), oob = squish) +
    scale_fill_gradientn(name = "Relative Interaction\nFrequency",
                         colours = c("#9E9E9E", "#ffffff", "#79248C"),
                         na.value = "white", 
                         limit = c(-1, 1),
                         # limit = c(0.95, 1.05),
                         oob = squish) +
    xlab("Individual") +
    ylab("Individual") +
    theme(axis.text = element_text(colour = "black", size = 6),
          axis.title = element_text(size = 7),
          axis.ticks = element_line(size = 0.3),
          aspect.ratio = 1,
          # Hide the legend (optional)
          legend.position = "none",
          legend.key.height = unit(0.38, "in"),
          panel.background = element_rect(size = 0.3, fill = NA),
          panel.grid = element_blank(),
          plot.title = element_blank()) +
    ggtitle(paste0("Group Size = ", groupsize))
  # return graph
  return(gg_avg_adj)
  # return(avg_g)
})


####################
# Single example group size plot
####################




####################
# Output example graph
####################
example_graph <- soc_networks[[14]][[1]]
example_thresh <- as.data.frame(thresh_data[[14]][[1]])
example_thresh$Id <- row.names(example_thresh)
example_thresh$ThreshRatio <- log(example_thresh$Thresh1 / example_thresh$Thresh2)
example_thresh$ThreshRatioBounded <- example_thresh$ThreshRatio
example_thresh$ThreshRatioBounded[example_thresh$ThreshRatioBounded < -5] <- -5
example_thresh$ThreshRatioBounded[example_thresh$ThreshRatioBounded > 5] <- 5
# Calculate values expected
reweight_graph <- example_graph
not_chosen <- 1 - (( 1 / (nrow(reweight_graph) - 1)) * p)
expected_random <-  1 - not_chosen^2
reweight_graph <- (reweight_graph - expected_random) / reweight_graph 
# Zero out interactions equal to or less than random
# example_graph[reweight_graph <= 0] <- 0
# Or take in those in top X percentile
percentiles <- quantile(example_graph, na.rm = TRUE)
fiftypercent <- percentiles[3]
seventyfivepercent <- percentiles[4]
example_graph[example_graph < fiftypercent] <- 0
diag(example_graph) <- 0
# Turn into graph object to get edgelist
g <- graph_from_adjacency_matrix(example_graph, mode = "undirected", weighted = TRUE)
edgelist <- get.edgelist(g)
edgelist <- as.data.frame(edgelist)
names(edgelist) <- c("Source", "Target")
edgelist$Weight <- E(g)$weight 
# Write
write.csv(edgelist, file = "output/Networks/ExampleNetworks/GroupSize70edgelist.csv", row.names = FALSE)
write.csv(example_thresh, file = "output/Networks/ExampleNetworks/GroupSize70nodelist.csv", row.names = FALSE)
