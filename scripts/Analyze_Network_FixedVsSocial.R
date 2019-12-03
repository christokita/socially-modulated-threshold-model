################################################################################
#
# Analyzing social network features comapred to activity and other factors
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

##########################################################
# DOL Entropy
##########################################################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

# Load data
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")
compiled_data$Model <- "Social"
entropy_data <- compiled_data
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0.05-Epsilon0-Beta1.1.Rdata")
compiled_data$Model <- "Fixed"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)

# Calculate mean and SE
entropy <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind),
            SE = sd(Dind) / sqrt(length(Dind)))

# Plot

gg_dol <- ggplot(data = entropy, aes(x = n, colour = Model, fill = Model)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                position = position_dodge(width = 1.5),
                size = 0.3) +
  geom_line(aes(y = Mean),
            position = position_dodge(width = 1.5),
            size = 0.4) +
  geom_point(aes(y = Mean),
             position = position_dodge(width = 1.5),
             size = 1, shape = 21) +
  theme_classic() +
  xlab(expression(paste("Group size (", italic(n), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20))) +
  scale_color_manual(name = "Threshold",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    values = c("#ffffff", "#4d4d4d")) +
  scale_linetype_manual(name = "Threshold",
                        values = c("dotted", "solid")) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = "none",
        legend.title = element_text(size = 7, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(5, "mm"),
        axis.ticks = element_line(size = 0.3),
        axis.line = element_line(size = 0.3),
        aspect.ratio = 1)
gg_dol

ggsave(gg_dol, filename = "output/SpecializationPlots/FixedVsSocialThresholds.svg", 
       height = 46, width = 46, units = "mm")


##########################################################
# Percentage of non-random interactions
##########################################################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")

p <- 1 #prob of interact
runs <- c("Sigma0.05-Epsilon0-Beta1.1", 
          "Sigma0-Epsilon0.1-Beta1.1")
run_names <- c("Fixed", "Social")

interaction_rates <- lapply(1:length(runs), function(run) {
  # Load social networks
  files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
  soc_networks <- list()
  for (file in 1:length(files)) {
    load(files[file])
    soc_networks[[file]] <- listed_data
  }
  # Load threshold matrices
  files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh/", runs[run], "/"), full.names = TRUE)
  thresh_data <- list()
  for (file in 1:length(files)) {
    load(files[file])
    thresh_data[[file]] <- listed_data
  }
  # Loop through individual graphs
  interaction_info <- lapply(1:length(soc_networks), function(i) {
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
      ratio <- order(thresh$ThreshBias)
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
    #  Calcualte 99% CI interval of interaction rate
    edgelist_sig <- all_edgelist %>% 
      group_by(From, To, Interaction) %>% 
      # filter(!is.na(Weight)) %>% 
      summarise(samp_mean = mean(Weight),
                samp_sd = sd(Weight),
                samples = length(Weight)) %>% 
      mutate(error = qt(0.975, df = samples-1) * samp_sd/sqrt(samples),
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
      summarise(Nonrandom = sum(DiffDirection!=0),
                HigherThanRandom = sum(DiffDirection == 1),
                LowerThanRandom = sum(DiffDirection == -1),
                TotalInteractions = n()) %>% 
      mutate(PercentNonRandom = Nonrandom / TotalInteractions,
             PercentHigher = HigherThanRandom / TotalInteractions,
             PercentLower = LowerThanRandom / TotalInteractions,
             n = dimensions[1],
             Model = run_names[run]) 
    # Return
    print(paste0(run_names[run], ": ", dimensions[1]))
    return(edgelist_sig)
  })
  # Bind and return
  interaction_info <- do.call("rbind", interaction_info)
  return(interaction_info)
})

# Bind
interaction_data <- do.call('rbind', interaction_rates)

# Graph
gg_interactions <- ggplot(interaction_data, aes(x = n, y = PercentNonRandom, 
                                                colour = Model, group = Model, fill = Model)) +
  geom_line(size = 0.4) +
  geom_point(size = 0.8, shape = 21) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(name = "Threshold type",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold type",
                    values = c("#ffffff", "#4d4d4d")) +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20))) +
  xlab(expression(paste("Group size ", italic(n)))) +
  ylab("% Non-random interactions") +
  # facet_grid(~Model) +
  theme_ctokita() +
  theme(aspect.ratio = 1,
        # strip.background = element_blank(),
        # strip.text = element_blank(),
        legend.position = "none")

gg_interactions
ggsave(gg_interactions, filename = "Output/Networks/NetworkMetrics/PercentNonRandomInteractions.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_interactions, filename = "Output/Networks/NetworkMetrics/PercentNonRandomInteractions.svg", 
       height = 45, width = 45, units = "mm")
#Facet plot
# ggsave(gg_interactions, filename = "Output/Networks/NetworkMetrics/PercentNonRandomInteractions_Facet.svg",
# width = 90, units = "mm")

##########################################################
# Modularity
##########################################################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")

p <- 1 #prob of interact
runs <- c("Sigma0.05-Epsilon0-Beta1.1",
          "Sigma0-Epsilon0.1-Beta1.1")
run_names <- c("Fixed", "Social")

modularity <- lapply(1:length(runs), function(run) {
  # Load social networks
  files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
  soc_networks <- list()
  for (file in 1:length(files)) {
    load(files[file])
    soc_networks[[file]] <- listed_data
  }
  # Load threshold matrices
  files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh/", runs[run], "/"), full.names = TRUE)
  thresh_data <- list()
  for (file in 1:length(files)) {
    load(files[file])
    thresh_data[[file]] <- listed_data
  }
  # Loop through individual graphs
  interaction_info <- lapply(1:length(soc_networks), function(i) {
    # Get graphs
    graphs <- soc_networks[[i]]
    replicates <- length(graphs)
    # For each each compute interaction matrix
    # Get graph and make adjacency matrix
    size_graph <- lapply(1:length(graphs), function(j) {
      # Format: set diagonal, rescale, and make adj matrix
      this_graph <- graphs[[j]]
      diag(this_graph) <- 0
      g <- graph_from_adjacency_matrix(this_graph, mode = "undirected", weighted = TRUE)
      g_clust <- cluster_fast_greedy(g, weights = E(g)$weight)
      # g_membership <- membership(g_clust)
      # mod <- modularity(g, membership = g_membership, weights = E(g)$weight)
      mod <- modularity(g_clust)
      clust_coeff <- transitivity(graph = g, type = "weighted", weights = E(g)$weight)
      # return
      replicate_row <- data.frame(n = nrow(this_graph),
                                  Modularity = mod,
                                  ClustCoeff =  mean(clust_coeff, na.rm = TRUE))
      return(replicate_row)
    })
    size_data <- do.call("rbind", size_graph)
  })
  # Bind and return
  interaction_info <- do.call("rbind", interaction_info)
  interaction_info$Model <- run_names[run]
  return(interaction_info)
})
# Bind
mod_data <- do.call("rbind", modularity)
mod_data <- mod_data %>% 
  group_by(Model, n) %>% 
  summarise(Modul_mean = mean(Modularity),
            Modul_SD = sd(Modularity),
            Modul_SE = sd(Modularity)/length(Modularity))

# Plot
gg_mod <- ggplot(mod_data, aes(x = n, y = Modul_mean, colour = Model, fill = Model)) +
  geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = Modul_mean - Modul_SD, ymax = Modul_mean + Modul_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  scale_color_manual(name = "Threshold",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    values = c("#ffffff", "#4d4d4d")) +
  scale_linetype_manual(name = "Threshold",
                        values = c("dotted", "solid")) +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20))) +
  scale_y_continuous(breaks = seq(0, 0.012, 0.002), limits = c(-0.0002, 0.012)) +
  xlab(expression(paste("Group size ", italic(n)))) +
  ylab("Modularity") +
  theme_ctokita() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        legend.key.height = unit(0.5, "line"))

gg_mod

ggsave(gg_mod, filename = "Output/Networks/NetworkMetrics/Modularity.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_mod, filename = "Output/Networks/NetworkMetrics/Modularity.svg", 
       height = 46, width = 46, units = "mm")

##########################################################
# Assortivity
##########################################################
runs <- c("Sigma0.05-Epsilon0-Beta1.1", 
          "Sigma0-Epsilon0.1-Beta1.1")
run_names <- c("Fixed", "Social")

###################
# Assortment coefficient from Newman 2003
###################
library(assortnet)

network_assort <- lapply(1:length(runs), function(run) {
  print(runs[run])
  # Load social networks
  files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
  soc_networks <- list()
  for (file in 1:length(files)) {
    load(files[file])
    soc_networks[[file]] <- listed_data
  }
  # Load threshold matrices
  files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh/", runs[run], "/"), full.names = TRUE)
  thresh_data <- list()
  for (file in 1:length(files)) {
    load(files[file])
    thresh_data[[file]] <- listed_data
  }
  # Loop through individual graphs
  interaction_info <- lapply(1:length(soc_networks), function(i) {
    print(i * 5)
    # Get graphs
    graphs <- soc_networks[[i]]
    replicates <- length(graphs)
    # For each each compute interaction matrix
    # Get graph and make adjacency matrix
    size_graph <- lapply(1:length(graphs), function(j) {
      # Get graph and calculate threshold differences
      this_graph <- graphs[[j]]
      diag(this_graph) <- 0
      thresh <- as.data.frame(thresh_data[[i]][j])
      thresh$ThreshBias <- thresh$Thresh1 - thresh$Thresh2
      # Calculate assortmnet
      assort <- assortment.continuous(graph = this_graph, vertex_values = thresh$ThreshBias, weighted = T)
      assort <- assort$r
      to_retun <- data.frame(n = nrow(this_graph), Assortativity = assort)
      # return
      return(to_retun)
    })
    #Calculate baseline probability of interaction
    size_graph <- do.call("rbind", size_graph)
  })
  # Bind and return
  interaction_info <- do.call("rbind", interaction_info)
  interaction_info$Model <- run_names[run]
  return(interaction_info)
})

# Bind
assort_data <- do.call('rbind', network_assort)
assort_data <- assort_data %>%
  group_by(Model, n) %>%
  summarise(Assort_mean = mean(Assortativity),
            Assort_SD = sd(Assortativity), 
            Assort_SE = sd(Assortativity)/length(Assortativity))

# Plot
gg_assort <- ggplot(data = assort_data, aes(x = n, y = Assort_mean,
                                                      colour = Model, group = Model, fill = Model)) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dotted") +
  geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = Assort_mean - Assort_SD, ymax = Assort_mean + Assort_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  scale_color_manual(name = "Threshold",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    values = c("#ffffff", "#4d4d4d")) +
  scale_linetype_manual(name = "Threshold",
                        values = c("dotted", "solid")) +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20))) +
  scale_y_continuous(breaks = seq(-0.25, 0.05, 0.05), limits = c(-0.26, 0.05)) +
  xlab(expression(paste("Group size ", italic(n)))) +
  ylab("Assortativity") +
  theme_ctokita() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        legend.key.height = unit(0.5, "line"))

gg_assort

ggsave(gg_assort, filename = "Output/Networks/NetworkMetrics/AssortmentCoeff.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_assort, filename = "Output/Networks/NetworkMetrics/AssortmentCoeff.svg", 
       height = 46, width = 46, units = "mm")



##########################################################
# Clustering
##########################################################
library(igraph)
library(tnet)

runs <- c("Sigma0.05-Epsilon0-Beta1.1",
          "Sigma0-Epsilon0.1-Beta1.1")
run_names <- c("Fixed", "Social")

###################
# Clustering coeff based on interactions above ranom  = 1 and all else  = 0
###################
# Note: tried normal clustering coeff before

network_clust <- lapply(1:length(runs), function(run) {
  print(runs[run])
  # Load social networks
  files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
  soc_networks <- list()
  for (file in 1:length(files)) {
    load(files[file])
    soc_networks[[file]] <- listed_data
  }
  # Load threshold matrices
  files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh/", runs[run], "/"), full.names = TRUE)
  thresh_data <- list()
  for (file in 1:length(files)) {
    load(files[file])
    thresh_data[[file]] <- listed_data
  }
  # Loop through individual graphs
  interaction_info <- lapply(1:length(soc_networks), function(i) {
    print(i * 5)
    # Get graphs
    graphs <- soc_networks[[i]]
    replicates <- length(graphs)
    # For each each compute interaction matrix
    # Get graph and make adjacency matrix
    size_graph <- lapply(1:length(graphs), function(j) {
      # Get graph and calculate baseline probability of interaction
      this_graph <- graphs[[j]]
      diag(this_graph) <- 0
      #Calculate baseline probability of interaction
      # dimensions <- dim(this_graph)
      # not_chosen <- 1 - (( 1 / (dimensions[1] - 1)))
      # expected_random <-  1 - not_chosen^2
      # # make all interactions above random = 1
      # this_graph[this_graph > expected_random] <- 1
      # this_graph[this_graph <= expected_random] <- 0
      # g <- graph_from_adjacency_matrix(this_graph)
      # # Calculate clustering
      # clust <- transitivity(graph = g, type = "global")
      # to_return <- data.frame(n = nrow(this_graph), ClusterCoeff = clust)
      
      # using Opsahl 2009 method
      g <- graph_from_adjacency_matrix(this_graph, weighted = T)
      e_list <- as.data.frame(get.edgelist(g))
      e_list <- e_list %>% 
        mutate(V1 = as.character(V1),
               V2 = as.character(V2)) %>% 
        mutate(V1 = as.numeric(gsub("v-", "", V1)),
               V2 = as.numeric(gsub("v-", "", V2)))
      e_list$weight <- E(g)$weight
      clust <- clustering_w(e_list)
      to_return <- data.frame(n = nrow(this_graph), ClusterCoeff = clust)
      # return
      return(to_return)
    })
    #Calculate baseline probability of interaction
    size_graph <- do.call("rbind", size_graph)
  })
  # Bind and return
  interaction_info <- do.call("rbind", interaction_info)
  interaction_info$Model <- run_names[run]
  return(interaction_info)
})

# Bind
clust_data <- do.call('rbind', network_clust)
clust_data <- clust_data %>%
  group_by(Model, n) %>%
  summarise(Clust_mean = mean(ClusterCoeff),
            Clust_SD = sd(ClusterCoeff),
            Clust_SE = sd(ClusterCoeff)/length(ClusterCoeff))

# Plot
gg_clust <- ggplot(data = clust_data, aes(x = n, y = Clust_mean,
                                            colour = Model, group = Model, fill = Model)) +
  geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = Clust_mean - Clust_SD, ymax = Clust_mean + Clust_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  scale_color_manual(name = "Threshold",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    values = c("#ffffff", "#4d4d4d")) +
  scale_linetype_manual(name = "Threshold",
                        values = c("dotted", "solid")) +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20))) +
  # scale_y_continuous(breaks = seq(-0.25, 0.05, 0.05), limits = c(-0.26, 0.05)) +
  xlab(expression(paste("Group size ", italic(n)))) +
  ylab("Clustering Coeff.") +
  theme_ctokita() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        legend.key.height = unit(0.5, "line"))

gg_clust





###################
# Weighted correlation (analogous to assortivity?)
###################
# library(wCorr)
# weighted_correlation <- lapply(1:length(runs), function(run) {
#   print(runs[run])
#   # Load social networks
#   files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
#   soc_networks <- list()
#   for (file in 1:length(files)) {
#     load(files[file])
#     soc_networks[[file]] <- listed_data
#   }
#   # Load threshold matrices
#   files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh/", runs[run], "/"), full.names = TRUE)
#   thresh_data <- list()
#   for (file in 1:length(files)) {
#     load(files[file])
#     thresh_data[[file]] <- listed_data
#   }
#   # Loop through individual graphs
#   interaction_info <- lapply(1:length(soc_networks), function(i) {
#     print(i*5)
#     # Get graphs
#     graphs <- soc_networks[[i]]
#     replicates <- length(graphs)
#     # For each each compute interaction matrix
#     # Get graph and make adjacency matrix
#     size_graph <- lapply(1:length(graphs), function(j) {
#       # Get graph and calculate threshold differences
#       this_graph <- graphs[[j]]
#       number_individuals <- dim(this_graph)[1]
#       diag(this_graph) <- 0
#       thresh <- as.data.frame(thresh_data[[i]][j])
#       thresh$ThreshBias <- thresh$Thresh1 - thresh$Thresh2 
#       # Multiply to get bias weighted by interaction frequency
#       social_interaction <- data.frame(ThreshBias = NULL, InteractBias = NULL, InteractWeight = NULL)
#       for (ind in 1:nrow(this_graph)) {
#         # Calculate bias and weighted interaction
#         thresh_bias <- rep(thresh$ThreshBias[ind], number_individuals)
#         interact_bias <- thresh$ThreshBias
#         interact_weight <- this_graph[ind, ]
#         interact_weight <- interact_weight[!is.na(interact_weight)]
#         # Bind to dataframe
#         to_bind <- data.frame(ThreshBias = thresh_bias, InteractBias = interact_bias, InteractWeight = interact_weight)
#         social_interaction <- rbind(social_interaction, to_bind)
#       }
#       # Filter and return
#       weighted_corr <- weightedCorr(x = social_interaction$ThreshBias,
#                                     y = social_interaction$InteractBias, 
#                                     method = "pearson", 
#                                     weights = social_interaction$InteractWeight)
#       to_return <- data.frame(n = number_individuals, WeightedCorr = weighted_corr)
#       # return
#       return(to_return)
#     })
#     #Calculate baseline probability of interaction
#     size_graph <- do.call("rbind", size_graph)
#   })
#   # Bind and return
#   interaction_info <- do.call("rbind", interaction_info)
#   interaction_info$Model <- run_names[run]
#   return(interaction_info)
# })
# 
# # Bind
# weightcorr_data <- do.call('rbind', weighted_correlation)
# weightcorr_data <- weightcorr_data %>% 
#   group_by(Model, n) %>% 
#   summarise(WeightedCorr_mean = mean(WeightedCorr),
#             WeightedCorr_SD = sd(WeightedCorr),
#             WeightedCorr_SE = sd(WeightedCorr)/length(WeightedCorr))
# 
# # Plot
# gg_weightedcorr <- ggplot(data = weightcorr_data, aes(x = n, y = WeightedCorr_mean, 
#                                                       colour = Model, group = Model, fill = Model)) +
#   geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dotted") +
#   geom_errorbar(aes(ymin = WeightedCorr_mean - WeightedCorr_SD, ymax = WeightedCorr_mean + WeightedCorr_SD),
#                 width = 0,
#                 size = 0.3) +
#   geom_line(size = 0.4) +
#   geom_point(size = 0.8, shape = 21) +
#   scale_color_manual(name = "Threshold",
#                      values = c("#878787", "#4d4d4d")) +
#   scale_fill_manual(name = "Threshold",
#                     values = c("#ffffff", "#4d4d4d")) +
#   scale_linetype_manual(name = "Threshold",
#                         values = c("dotted", "solid")) +
#   scale_x_continuous(breaks = seq(0, 100, 20)) +
#   scale_y_continuous(breaks = seq(-0.25, 0.05, 0.05), limits = c(-0.26, 0.05)) +
#   xlab(expression(paste("Group size ", italic(n)))) +
#   ylab("Assortativity") +
#   theme_ctokita() +
#   theme(aspect.ratio = 1,
#         legend.position = "none",
#         legend.key.height = unit(0.5, "line"))
# gg_weightedcorr
# 
# ggsave(gg_weightedcorr, filename = "Output/Networks/NetworkMetrics/WeightedCorrelation.png", 
#        height = 45, width = 45, units = "mm", dpi = 400)
# ggsave(gg_weightedcorr, filename = "Output/Networks/NetworkMetrics/WeightedCorrelation.svg", 
#        height = 46, width = 46, units = "mm")

