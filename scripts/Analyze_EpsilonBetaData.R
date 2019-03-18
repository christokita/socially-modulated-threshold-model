################################################################################
#
# Analyze beta-epsilon sweep data
#
################################################################################

##############################  Entropy plots  ############################## 
rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

####################
# Load data
####################
load("output/ParameterSpace/EpsilonBetaSweep-n60.Rdata")

# Filter to parameter combos of interest
epsilons_of_interest <- c(0, 0.1, 0.3, 0.5)
betas_of_interest <- c(1, 1.05, 1.1, 1.2)

betas <- entropy %>% 
  filter(epsilon %in% epsilons_of_interest) %>% 
  mutate(Mean = Dind_mean,
         SD = Dind_SD,
         epsilon = as.factor(epsilon))

epsilons <- entropy %>% 
  filter(beta %in% betas_of_interest) %>% 
  mutate(Mean = Dind_mean,
         SD = Dind_SD,
         beta = as.factor(beta))

####################
# Plot betas
####################
# Betas
pal <- brewer.pal(5, "Greens")[2:5]

gg_entropy_betas <- ggplot(data = betas, aes(x = beta, colour = epsilon)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_color_manual(values = pal, 
                     labels = c("0.0", "0.1", "0.3", "0.5"),
                     name = expression("Social\ninfluence "(epsilon))) +
  theme_ctokita() +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = "none",
        legend.title = element_text(size = 6, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_entropy_betas

ggsave(gg_entropy_betas, file = "output/SpecializationPlots/Beta-Eps-n60-DiffEpsValues.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_entropy_betas, file = "output/SpecializationPlots/Beta-Eps-n60-DiffEpsValues.svg", 
       heigh = 45, width = 45, units = "mm")

# resize
gg_entropy_betas <- gg_entropy_betas + 
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme(aspect.ratio = NULL,
        axis.title.y = element_blank())

ggsave(gg_entropy_betas, file = "output/SpecializationPlots/Beta-Eps-n60-DiffEpsValues.svg", 
       width = 45, height = 25, units = "mm")

####################
# Plot epsilons
####################
# Epsilon
pal <- brewer.pal(5, "Greens")[2:5]

gg_entropy_eps <- ggplot(data = epsilons, aes(x = epsilon, colour = beta)) +
  geom_errorbar(aes(ymin = ifelse((Mean - SD) > 0 , Mean - SD, 0), ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Social influence (", italic(epsilon), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(0, 0.6, 0.1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = pal, 
                     labels = as.character(betas_of_interest),
                     name = expression("Interaction\nbias "(beta))) +
  theme_ctokita() +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = "right",
        legend.title = element_text(size = 6, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.ticks = element_line(size = 0.3),
        axis.line = element_line(size = 0.3),
        aspect.ratio = 1)
gg_entropy_eps

ggsave(gg_entropy_eps, file = "output/SpecializationPlots/Beta-Eps-n60-DiffBetaValues.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_entropy_eps, file = "output/SpecializationPlots/Beta-Eps-n60-DiffBetaValues.svg", 
       width = 45, height = 45, units = "mm")

# resize
gg_entropy_eps <- gg_entropy_eps + 
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme(aspect.ratio = NULL,
        axis.title.y = element_blank())

ggsave(gg_entropy_eps, file = "output/SpecializationPlots/Beta-Eps-n60-DiffBetaValues.svg", 
       width = 45, height = 25, units = "mm")


##############################  Network property plots  ############################## 

##########################################################
# Modularity
##########################################################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")

p <- 1 #prob of interact
runs <- c("Sigma0-Epsilon0.1_BetaSweep")
run_names <- c("Social")

modularity <- lapply(1:length(runs), function(run) {
  # Load social networks
  files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
  betas <- gsub(".*/([\\.0-9]+).Rdata", "\\1", x = files, perl = T)
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
      replicate_row <- data.frame(beta = as.numeric(betas[i]),
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
  group_by(Model, beta) %>% 
  summarise(Modul_mean = mean(Modularity),
            Modul_SD = sd(Modularity),
            Modul_SE = sd(Modularity)/length(Modularity))

# Plot
gg_mod <- ggplot(mod_data, aes(x = beta, y = Modul_mean, colour = Model, fill = Model)) +
  # geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = Modul_mean - Modul_SD, ymax = Modul_mean + Modul_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  scale_color_manual(name = "Threshold",
                     # values = c("#878787", "#4d4d4d")) +
                     values = c("#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    # values = c("#ffffff", "#4d4d4d")) +
                    values = c("#4d4d4d")) +
  scale_x_continuous(breaks = seq(1, 1.5, 0.05)) +
  scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(-0.0002, 0.03)) +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab("Modularity") +
  theme_ctokita() +
  theme(legend.position = "none",
        legend.key.height = unit(0.5, "line"))

gg_mod

# ggsave(gg_mod, filename = "Output/Networks/NetworkMetrics/Modularity_betasweep.svg", 
#        height = 23, width = 46, units = "mm")

#########################################################
# Assortivity
##########################################################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")

p <- 1 #prob of interact
runs <- c("Sigma0-Epsilon0.1_BetaSweep")
run_names <- c("Social")

###################
# Assortment coefficient from Newman 2003
###################
library(assortnet)
library(gridExtra)

network_assort <- lapply(1:length(runs), function(run) {
  print(runs[run])
  # Load social networks
  files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
  betas <- gsub(".*/([\\.0-9]+).Rdata", "\\1", x = files, perl = T)
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
      # Get graph and calculate threshold differences
      this_graph <- graphs[[j]]
      diag(this_graph) <- 0
      thresh <- as.data.frame(thresh_data[[i]][j])
      thresh$ThreshBias <- thresh$Thresh1 - thresh$Thresh2
      # Calculate assortmnet
      assort <- assortment.continuous(graph = this_graph, vertex_values = thresh$ThreshBias, weighted = T)
      assort <- assort$r
      to_retun <- data.frame(beta = as.numeric(betas[i]), Assortativity = assort)
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
  group_by(Model, beta) %>%
  summarise(Assort_mean = mean(Assortativity),
            Assort_SD = sd(Assortativity), 
            Assort_SE = sd(Assortativity)/length(Assortativity))

# Plot
gg_assort <- ggplot(data = assort_data, aes(x = beta, y = Assort_mean,
                                            colour = Model, group = Model, fill = Model)) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dotted") +
  geom_errorbar(aes(ymin = Assort_mean - Assort_SD, ymax = Assort_mean + Assort_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  scale_color_manual(name = "Threshold",
                     # values = c("#878787", "#4d4d4d")) +
                     values = c("#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    # values = c("#ffffff", "#4d4d4d")) +
                    values = c("#4d4d4d")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(-0.04, 0.1, 0.02), limits = c(-0.02, 0.06)) +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab("Assortativity") +
  theme_ctokita() +
  theme(legend.position = "none")

gg_assort

# ggsave(gg_assort, filename = "Output/Networks/NetworkMetrics/Assortativity_betasweep.svg", 
#        height = 23, width = 49, units = "mm")

# Together
gg_net_metrics <- grid.arrange(gg_mod, gg_assort, nrow = 2)
ggsave(gg_net_metrics, filename = "Output/Networks/NetworkMetrics/NetworkMetrics_betasweep.svg", 
       height = 45, width = 49, units = "mm")
