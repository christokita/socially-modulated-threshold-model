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
load("output/Rdata/Sigma0.01-FIXED-Bias1.1.Rdata")

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


####################
# Compare social network features
####################

###### Network "dispersion" (? - standard deviation over mean degree) ###### 
soc_dispersion <- lapply(soc_graphs, function(graph) {
  # Group size
  n <- nrow(graph)
  # Degrees
  degrees <- rowSums(graph)
  deg_mean <- mean(degrees)
  deg_sd <- sd(degrees)
  # compile and return
  to_return <- data.frame(n = n, DegreeMean = deg_mean, DegreeSD = deg_sd)
})
dispersion <- do.call("rbind", soc_dispersion)
dispersion <- dispersion %>% 
  mutate(Dispersion = DegreeSD / DegreeMean)
dispersionSummary <- dispersion %>% 
  mutate(Model = "Social") %>% 
  group_by(n, Model) %>% 
  summarise(Dispersion = mean(Dispersion),
            DegreeMean = mean(DegreeMean),
            DegreeSD = mean(DegreeSD))

fixed_dispersion <- lapply(fix_graphs, function(graph) {
  # Group size
  n <- nrow(graph)
  # Degrees
  degrees <- rowSums(graph)
  deg_mean <- mean(degrees)
  deg_sd <- sd(degrees)
  # compile and return
  to_return <- data.frame(n = n, DegreeMean = deg_mean, DegreeSD = deg_sd)
})
fixed_dispersion <- do.call("rbind", fixed_dispersion)
fixed_dispersion <- fixed_dispersion %>% 
  mutate(Dispersion = DegreeSD / DegreeMean)
fixed_dispersionSummary <- fixed_dispersion %>% 
  mutate(Model = "Fixed") %>% 
  group_by(n, Model) %>% 
  summarise(Dispersion = mean(Dispersion),
            DegreeMean = mean(DegreeMean),
            DegreeSD = mean(DegreeSD))

dispersionSummary <- rbind(dispersionSummary, fixed_dispersionSummary)
rm(fixed_dispersion, fixed_dispersionSummary)

gg_dispersion <- ggplot(data = dispersionSummary, aes(x = n, y = Dispersion, group = Model, color = Model)) +
  geom_line(aes(linetype = Model)) +
  geom_point(size = 2) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("black", "mediumseagreen")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black")) 
gg_dispersion


###### Network homophily ###### 
weighted_ratios <- lapply(1:length(soc_graphs), function(i) {
  # Social
  # Get graph and thresh matrix for simulation
  graph <- soc_graphs[[i]]
  thresh <- as.data.frame(soc_threshMat[[i]])
  # Calculate thresh ratio
  threshRatio <- log(thresh$Thresh1 / thresh$Thresh2)
  # Calculate weighted neighbor sum for each individual
  weighted_sum <- graph %*% threshRatio
  # Calculate weighted neighbor sum for each individual
  strength <- rowSums(graph)
  disparity <- rowSums((graph / strength)^2)
  # Construct dataframe to return
  to_return <- data.frame(Id = row.names(weighted_sum), 
                          n = length(weighted_sum),
                          ThreshRatio = threshRatio,
                          WeightNeighbor = weighted_sum,
                          Disparity = disparity,
                          Model = "Social")
  row.names(to_return) <- NULL
  # Fixed
  # Get graph and thresh matrix for simulation
  graph <- fix_graphs[[i]]
  thresh <- as.data.frame(fix_threshMat[[i]])
  # Calculate thresh ratio
  threshRatio <- log(thresh$Thresh1 / thresh$Thresh2)
  # Calculate weighted neighbor sum for each individual
  weighted_sum <- graph %*% threshRatio
  # Calculate weighted neighbor sum for each individual
  strength <- rowSums(graph)
  disparity <- rowSums((graph / strength)^2)
  # Construct dataframe to return
  to_return2 <- data.frame(Id = row.names(weighted_sum), 
                          n = length(weighted_sum),
                          ThreshRatio = threshRatio,
                          WeightNeighbor = weighted_sum,
                          Disparity = disparity,
                          Model = "Fixed")
  row.names(to_return) <- NULL
  # Return
  to_return$Model <- as.character(to_return$Model)
  to_return <- rbind(to_return, to_return2)
  return(to_return)
})
weighted_ratios <- do.call("rbind", weighted_ratios)

test <- weighted_ratios %>% filter(Model == "Social")

gg_weighted_ratios <- ggplot(data = test, aes(x = ThreshRatio, y = WeightNeighbor)) +
  geom_point() +
  theme_classic(base_size = 10) +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"))  +
  facet_wrap(~n, scales = "free_y")
gg_weighted_ratios

weighted_summary <- weighted_ratios %>% 
  group_by(n, Model) %>% 
  summarise(CorrelationThresholds = cor(ThreshRatio, WeightNeighbor))

gg_weighted_corr <- ggplot(data = weighted_summary, aes(x = n, y = CorrelationThresholds, group = Model, color = Model)) +
  geom_line(aes(linetype = Model)) +
  geom_point(size = 2) +
  theme_classic(base_size = 10) +
  scale_color_manual(values = c("black", "mediumseagreen")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black")) 
gg_weighted_corr


###### Network assortivity ######
# From Newmann (2003)
assortivity <- lapply(1:length(soc_graphs), function(i) {
  # Social
  # Get graph and thresh matrix for simulation
  graph <- soc_graphs[[i]]
  thresh <- as.data.frame(soc_threshMat[[i]])
  # Calculate thresh ratio
  threshRatio <- log(thresh$Thresh1 / thresh$Thresh2)
  # Normalize graph
  norm_graph <- graph / sum(graph)
  # Calculate
  sigma_a <- sd(rowSums(norm_graph))
  pearson <- lapply(1:length(threshRatio), function(n) {
    x <- threshRatio[n]
    to_sum <- c()
    for (m in n:length(threshRatio)) {
      y <- threshRatio[m]
      e_xy <- norm_graph[n, m]  #symmetric
      a_x <- sum(norm_graph[n, ])
      value <- x * y * (e_xy - 2 * a_x)
      to_sum <- c(to_sum, value)
    }
    sum_value <- sum(to_sum)
    return(sum_value)
  })
  pearson <- sum(unlist(pearson)) / sigma_a / sigma_a
  
})
assortivity <- do.call("rbind", assortivity)




###### Network disparity ###### 
dis_test <- weighted_ratios %>% filter(Model == "Fixed")

gg_disparity <- ggplot(data = dis_test, aes(x = ThreshRatio, y = Disparity)) +
  geom_point(size = 0.2) +
  theme_classic(base_size = 10) +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"))  +
  facet_wrap(~n, scales = "free")
gg_disparity


###### Interaction Matrix Plot ###### 

type_groups_graphs <- soc_groups_graphs
type_groups_threshMat <- soc_groups_threshMat

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
    this_graph <- scale(this_graph)
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
                         colours = brewer.pal(5,"BuPu"), 
                         na.value = "white", 
                         limit = c(-0.5, 0.5),
                         oob = squish) +
    theme(axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      aspect.ratio = 1,
      # Hide the legend (optional)
      legend.position = "right",
      legend.key.height = units(2, "in"),
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
png(filename = paste0("output/Networks/FIGURE_IntMat_Single", filename, ".png"), width = 2, height = 2, units = "in", res = 600)
multiplot(plotlist = interaction_graphs[c(7)], layout = matrix(c(seq(1:length(interaction_graphs[c(1, 2, 5, 8)]))), nrow = 1, byrow=TRUE))
dev.off()


###### Distribution of degree ###### 
degree_dist <- lapply(1:length(soc_groups_graphs), function(i) {
  # Get graphs
  graphs <- soc_groups_graphs[[i]]
  replicates <- length(graphs)
  
})
