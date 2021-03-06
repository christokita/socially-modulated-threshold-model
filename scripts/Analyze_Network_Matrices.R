################################################################################
#
# Analyzing average social networks of simulations
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)
library(viridis)

p <- 1 #prob of interact
run <- "Sigma0.05-Epsilon0-Beta1.1"

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
    # Calculate thresh bias
    thresh <- as.data.frame(thresh_data[[i]][j])
    thresh$ThreshBias <- thresh$Thresh1 - thresh$Thresh2
    bias <- order(thresh$ThreshBias) # sort from most negative (bias task 1) to most positive (bias task 2)
    this_graph <- this_graph[bias, bias]
    colnames(this_graph) <- paste0("i-", 1:dimensions[1])
    rownames(this_graph) <- paste0("i-", 1:dimensions[1])
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
  if (groupsize < 20) {
    breaks <- c(1, seq(5, length(unique(plot_data$to)), 5))
  } else if(groupsize < 35) {
    breaks <- c(1, seq(10, length(unique(plot_data$to)), 15))
  } else {
    breaks <- c(1, seq(20, length(unique(plot_data$to)), 20))
  }
  # Color palette
  # pal <- c('#525252','#6c6c6c','#878787','#a4a4a4','#c2c2c2','#e0e0e0', 
  #          '#ffffff',
  #          '#e2d7eb','#c7b1d7','#ad8ac1','#9763aa','#823b8f','#6e016b')
  pal <- c('#616161','#939393','#c7c7c7',
           '#ffffff',
           '#e1b9d1','#c074a5','#9b287b')
  # Plot
  gg_avg_adj <- ggplot(plot_data, aes(x = from, y = to, fill = weight, color = weight)) +
    geom_tile() +
    theme_bw() +
    # Because we need the x and y axis to display every node,
    # not just the nodes that have connections to each other,
    # make sure that ggplot does not drop unused factor levels
    scale_x_discrete(drop = FALSE, expand = c(0, 0), 
                     position = "top",
                     breaks = levels(plot_data$to)[breaks],
                     labels = rep("", length(levels(plot_data$to)[breaks]))) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0), 
                     limits = rev(levels(plot_data$to)),
                     breaks = levels(plot_data$to)[breaks],
                     labels = rep("", length(levels(plot_data$to)[breaks]))) +
    scale_fill_gradientn(name = "Relative\ninteraction\nfrequency",
                         colours = pal,
                         na.value = "white",
                         limit = c(-0.04, 0.04),
                         oob = squish) +
    scale_color_gradientn(name = "Relative\ninteraction\nfrequency",
                         colours =  pal,
                         na.value = "white",
                         limit = c(-0.04, 0.04),
                         oob = squish) +
    xlab("Individual") +
    ylab("Individual") +
    theme(aspect.ratio = 1,
          # Hide the legend (optional)
          legend.key.width = unit(3, "mm"),
          legend.key.height = unit(6, "mm"),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          panel.border = element_rect(size = 0.3, fill = NA, colour = "black"),
          plot.title = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          # # Minimal settings
          # axis.ticks = element_blank(),
          # Detailed Settings
          axis.ticks = element_line(size = 0.3, colour = "black"),
          legend.position = "none") +
    ggtitle(paste0("Group Size = ", groupsize))
  # return graph
  return(gg_avg_adj)
  # return(avg_g)
})

# Save (only for minimal plots)
plots <- seq(4, 7, 1)
for (plot in plots) {
  gg_inter <- interaction_graphs[[plot]]
  ggsave(gg_inter, 
         filename = paste0("output/Networks/RawPlots/", run, "_n", plot*5 ,".svg"), 
         width = 20, 
         height = 20,
         units = "mm")
}

# save specific plot (only for fully detailed plot)
gg_inter <- interaction_graphs[80/5]
gg_inter
ggsave("output/Networks/RawPlots/GroupSize80.svg", width = 38, height = 38, units = "mm")


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
    thresh$ThreshBias <- thresh$Thresh1 - thresh$Thresh2
    bias <- order(thresh$ThreshBias)
    # Create order by threshold bias
    this_graph <- this_graph[bias, bias]
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
  #  Calcualte 99% CI interval of interaction rate
  edgelist_sig <- all_edgelist %>% 
    group_by(From, To, Interaction) %>% 
    # filter(!is.na(Weight)) %>% 
    summarise(samp_mean = mean(Weight),
              samp_sd = sd(Weight),
              samples = length(Weight)) %>% 
    mutate(error = qt(0.975,df = samples-1) * samp_sd/sqrt(samples),
           CI_low = samp_mean - error,
           CI_high = samp_mean + error) %>% 
    mutate(Lower_check = CI_low > expected_random,
           Higher_check = CI_high > expected_random) 
  # Determine if it is different than random
  edgelist_sig <- as.data.frame(edgelist_sig)
  edgelist_sig$DiffDirection <- 0
  edgelist_sig$DiffDirection[edgelist_sig$Lower_check & edgelist_sig$Higher_check] <- 1
  edgelist_sig$DiffDirection[edgelist_sig$Lower_check==FALSE & edgelist_sig$Higher_check==FALSE] <- -1
  # Make graph
  edgelist_sig <- edgelist_sig %>% 
    select(From, To, Interaction, DiffDirection) %>% 
    # unique(.) %>% 
    mutate(From = as.numeric(as.character(From)),
           To = as.numeric(as.character(To))) %>% 
    arrange(From, To)
  node_list <- unique( paste0("i-", edgelist_sig$To))
  plot_data <- edgelist_sig %>% 
    mutate(from = paste0("i-", From),
           to = paste0("i-", To)) %>% 
    mutate(to = factor(to, levels = rev(node_list)),
           from = factor(from, levels = node_list),
           weight = DiffDirection)
  # Get info for plot
  groupsize <- length(node_list)
  if (groupsize < 30) {
    breaks <- c(1, seq(5, groupsize, 5))
  } else if(groupsize < 55) {
    breaks <- c(1, seq(10, groupsize, 10))
  } else {
    breaks <- c(1, seq(20, groupsize, 20))
  }
  # Plot
  gg_avg_adj <- ggplot(plot_data, aes(x = from, y = to, fill = weight, colour = weight)) +
    geom_tile() +
    theme_bw() +
    # Because we need the x and y axis to display every node,
    # not just the nodes that have connections to each other,
    # make sure that ggplot does not drop unused factor levels
    scale_x_discrete(drop = FALSE, expand = c(0, 0), 
                     position = "top",
                     breaks = levels(plot_data$from)[breaks]) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0), 
                     limits = levels(plot_data$to),
                     breaks = levels(plot_data$from)[breaks]) +
    scale_fill_gradientn(name = "Relative Interaction\nFrequency",
                         colours = c("#616161", "#ffffff", "#9b287b"),
                         na.value = "white",
                         limit = c(-1, 1),
                         oob = squish) +
    scale_color_gradientn(name = "Relative Interaction\nFrequency",
                         colours = c("#616161", "#ffffff", "#9b287b"),
                         na.value = "white",
                         limit = c(-1, 1),
                         oob = squish) +
    xlab("Individual") +
    ylab("Individual") +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          # axis.ticks = element_line(size = 0.3, colour = "black"),
          axis.ticks = element_blank(),
          # Hide the legend (optional)
          legend.position = "none",
          legend.key.width = unit(3, "mm"),
          legend.key.height = unit(6, "mm"),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          panel.border = element_rect(size = 0.3, fill = NA, colour = "black"),
          plot.title = element_blank(),
          aspect.ratio = 1,
          panel.grid = element_blank()) +
    ggtitle(paste0("Group Size = ", groupsize))
  # return graph
  print(groupsize)
  return(gg_avg_adj)
})

# Save single plot
gg_simp <- simple_graphs[80/5]
gg_simp
ggsave("output/Networks/RawPlots/SimpleAdjPlot_80.svg", width = 38, height = 38, units = "mm")

# Save all for gif
new_dir <- paste0("output/Networks/RawPlots/", run, "/")
dir.create(new_dir, showWarnings = FALSE)
for (i in 1:length(simple_graphs)) {
  gg_plot <- simple_graphs[[i]]
  plot_name <- i*5
  ggsave(gg_plot, filename = paste0(new_dir, plot_name, ".svg"), width = 15, height = 15, units = "mm")
}


