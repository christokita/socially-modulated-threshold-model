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
# Plot betas
####################
# Load data and summarise
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep.Rdata")

betas <- compiled_data %>% 
  group_by(beta) %>% 
  summarise(Mean = mean(Dind),
             SD = sd(Dind))

# Betas
gg_entropy_betas <- ggplot(data = betas, aes(x = beta)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3, 
                color = "#4d4d4d") +
  geom_point(aes(y = Mean),
             size = 0.8, 
             color = "#4d4d4d") +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme_ctokita() +
  theme(axis.title.y = element_text(vjust = -1.5))
gg_entropy_betas

ggsave(gg_entropy_betas, file = "output/SpecializationPlots/BetaSweep-n80-ForNetworkPlot.png", 
       height = 23, width = 45, units = "mm")


####################
# Plot epsilons
####################
# Load data and summarise
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Beta1.1_EpsSweep.Rdata")

epsilons <- compiled_data %>% 
  group_by(epsilon) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

gg_entropy_eps <- ggplot(data = epsilons, aes(x = epsilon)) +
  geom_errorbar(aes(ymin = ifelse((Mean - SD) > 0 , Mean - SD, 0), ymax = Mean + SD),
                width = 0,
                size = 0.3, 
                color = "#4d4d4d") +
  geom_point(aes(y = Mean),
             size = 0.8,
             color = "#4d4d4d") +
  theme_classic() +
  xlab(expression(paste("Social influence (", italic(epsilon), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(0, 0.6, 0.1)) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.5)) +
  theme_ctokita() +
  theme(axis.title.y = element_text(vjust = -1.5))
gg_entropy_eps

ggsave(gg_entropy_eps, file = "output/SpecializationPlots/EpsSweep-n80-ForNetworkPlot.png", 
       height = 23, width = 45, units = "mm")

####################
# Plot group size
####################
# Load data and summarise
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")

groupsizes <- compiled_data %>% 
  group_by(n) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

gg_entropy_gs <- ggplot(data = groupsizes, aes(x = n)) +
  geom_errorbar(aes(ymin = ifelse((Mean - SD) > 0 , Mean - SD, 0), ymax = Mean + SD),
                width = 0,
                size = 0.3, 
                color = "#4d4d4d") +
  geom_point(aes(y = Mean),
             size = 0.8,
             color = "#4d4d4d") +
  theme_classic() +
  xlab(expression(paste("Group size (", italic(n), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.5)) +
  theme_ctokita() +
  theme(axis.title.y = element_text(vjust = -1.5))
gg_entropy_gs

ggsave(gg_entropy_gs, file = "output/SpecializationPlots/GroupSizeSweep-Beta1.1Eps0.1-ForNetworkPlot.svg", 
       height = 23, width = 45, units = "mm")


####################
# Plot together
####################
library(gridExtra)
gg_entropy_betaeps <- grid.arrange(gg_entropy_betas, gg_entropy_eps, nrow = 1)
ggsave(gg_entropy_betaeps, filename = "Output/SpecializationPlots/DOL_betaepssweepForNetFig.svg", 
       height = 23, width = 93.75, units = "mm")

gg_entropy_togther <- grid.arrange(gg_entropy_betas, gg_entropy_eps, gg_entropy_gs, nrow = 1)
ggsave(gg_entropy_togther, filename = "Output/SpecializationPlots/DOL_parametersweepForNetFig.svg", 
       height = 23.1, width = 1.515*93.75, units = "mm")



##############################  Network property plots  ############################## 

##########################################################
# Modularity
##########################################################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")

p <- 1 #prob of interact
runs <- c("Sigma0-Epsilon0.1_BetaSweep",
          "Sigma0-Beta1.1_EpsSweep",
          "Sigma0-Epsilon0.1-Beta1.1")
run_names <- c("Beta", "Epsilon", "GroupSize")

modularity <- lapply(1:length(runs), function(run) {
  # Load social networks
  files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
  parameter_value <- gsub(".*/([\\.0-9]+).Rdata", "\\1", x = files, perl = T)
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
      replicate_row <- data.frame(parameter_value = as.numeric(parameter_value[i]),
                                  Modularity = mod,
                                  ClustCoeff =  mean(clust_coeff, na.rm = TRUE))
      return(replicate_row)
    })
    size_data <- do.call("rbind", size_graph)
  })
  # Bind and return
  interaction_info <- do.call("rbind", interaction_info)
  interaction_info$parameter <- run_names[run]
  return(interaction_info)
})
# Bind
mod_data <- do.call("rbind", modularity)
mod_data <- mod_data %>% 
  group_by(parameter, parameter_value) %>% 
  summarise(Modul_mean = mean(Modularity),
            Modul_SD = sd(Modularity),
            Modul_SE = sd(Modularity)/length(Modularity))

# Plot
mod_data_beta <- mod_data %>% 
  filter(parameter == "Beta")
gg_mod_beta <- ggplot(mod_data_beta, aes(x = parameter_value, y = Modul_mean, colour = parameter, fill = parameter)) +
  # geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = ifelse(Modul_mean - Modul_SD < 0, 0, Modul_mean - Modul_SD), ymax = Modul_mean + Modul_SD),
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
  scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(-0.0002, 0.031)) +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab("Modularity") +
  theme_ctokita() +
  theme(legend.position = "none",
        legend.key.height = unit(0.5, "line"))
gg_mod_beta

mod_data_eps <- mod_data %>% 
  filter(parameter == "Epsilon")
gg_mod_eps <- ggplot(mod_data_eps, aes(x = parameter_value, y = Modul_mean, colour = parameter, fill = parameter)) +
  # geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = ifelse(Modul_mean - Modul_SD < 0, 0, Modul_mean - Modul_SD), ymax = Modul_mean + Modul_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  scale_color_manual(name = "Threshold",
                     # values = c("#878787", "#4d4d4d")) +
                     values = c("#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    # values = c("#ffffff", "#4d4d4d")) +
                    values = c("#4d4d4d")) +
  scale_x_continuous(breaks = seq(0, 0.6, 0.1)) +
  scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(-0.0002, 0.031)) +
  xlab(expression(paste("Social influence (", italic(epsilon), ")"))) +
  ylab("Modularity") +
  theme_ctokita() +
  theme(legend.position = "none",
        legend.key.height = unit(0.5, "line"))
gg_mod_eps


mod_data_gs <- mod_data %>% 
  filter(parameter == "GroupSize")
gg_mod_gs <- ggplot(mod_data_gs, aes(x = parameter_value, y = Modul_mean, colour = parameter, fill = parameter)) +
  # geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = ifelse(Modul_mean - Modul_SD < 0, 0, Modul_mean - Modul_SD), ymax = Modul_mean + Modul_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  scale_color_manual(name = "Threshold",
                     # values = c("#878787", "#4d4d4d")) +
                     values = c("#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    # values = c("#ffffff", "#4d4d4d")) +
                    values = c("#4d4d4d")) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(-0.0002, 0.031)) +
  xlab(expression(paste("Group size (", italic(n), ")"))) +
  ylab("Modularity") +
  theme_ctokita() +
  theme(legend.position = "none",
        legend.key.height = unit(0.5, "line"))
gg_mod_gs



# ggsave(gg_mod, filename = "Output/Networks/NetworkMetrics/Modularity_betasweep.svg", 
#        height = 23, width = 46, units = "mm")

#########################################################
# Assortivity
##########################################################
source("scripts/util/__Util__MASTER.R")

p <- 1 #prob of interact
runs <- c("Sigma0-Epsilon0.1_BetaSweep",
          "Sigma0-Beta1.1_EpsSweep",
          "Sigma0-Epsilon0.1-Beta1.1")
run_names <- c("Beta", "Epsilon", "GroupSize")

###################
# Assortment coefficient from Newman 2003
###################
library(assortnet)
library(gridExtra)

network_assort <- lapply(1:length(runs), function(run) {
  print(runs[run])
  # Load social networks
  files <- list.files(paste0("output/Rdata/_ProcessedData/Graphs/", runs[run], "/"), full.names = TRUE)
  parameter_value <- gsub(".*/([\\.0-9]+).Rdata", "\\1", x = files, perl = T)
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
      to_retun <- data.frame(parameter_value = as.numeric(parameter_value[i]), Assortativity = assort)
      # return
      return(to_retun)
    })
    #Calculate baseline probability of interaction
    size_graph <- do.call("rbind", size_graph)
  })
  # Bind and return
  interaction_info <- do.call("rbind", interaction_info)
  interaction_info$parameter <- run_names[run]
  return(interaction_info)
})

# Bind
assort_data <- do.call('rbind', network_assort)
assort_data <- assort_data %>%
  mutate(Assortativity = ifelse(Assortativity == -1, NA, Assortativity)) %>% #remove huge outlier
  group_by(parameter, parameter_value) %>%
  summarise(Assort_mean = mean(Assortativity, na.rm = T),
            Assort_SD = sd(Assortativity, na.rm = T), 
            Assort_SE = sd(Assortativity, na.rm = T)/length(Assortativity))

# Plot
assort_data_beta <- assort_data %>% 
  filter(parameter == "Beta")
gg_assort_beta <- ggplot(data = assort_data_beta, aes(x = parameter_value, y = Assort_mean,
                                            colour = parameter, group = parameter, fill = parameter)) +
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
  scale_y_continuous(breaks = seq(-0.04, 0.1, 0.02), limits = c(-0.02, 0.062)) +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab("Assortativity") +
  theme_ctokita() +
  theme(legend.position = "none",
        axis.title.y = element_text(vjust = -1.5))
gg_assort_beta

assort_data_eps <- assort_data %>% 
  filter(parameter == "Epsilon")
gg_assort_eps <- ggplot(data = assort_data_eps, aes(x = parameter_value, y = Assort_mean,
                                                      colour = parameter, group = parameter, fill = parameter)) +
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
  scale_x_continuous(breaks = seq(0, 0.6, 0.1)) +
  scale_y_continuous(breaks = seq(-0.04, 0.1, 0.02), limits = c(-0.02, 0.062)) +
  xlab(expression(paste("Social influence (", italic(epsilon), ")"))) +
  ylab("Assortativity") +
  theme_ctokita() +
  theme(legend.position = "none",
        axis.title.y = element_text(vjust = -1.5))
gg_assort_eps

assort_data_gs <- assort_data %>% 
  filter(parameter == "GroupSize")
gg_assort_gs <- ggplot(data = assort_data_gs, aes(x = parameter_value, y = Assort_mean,
                                                    colour = parameter, group = parameter, fill = parameter)) +
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
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  # scale_y_continuous(breaks = seq(-0.04, 0.1, 0.02), limits = c(-0.02, 0.062)) +
  xlab(expression(paste("Group size (", italic(n), ")"))) +
  ylab("Assortativity") +
  theme_ctokita() +
  theme(legend.position = "none",
        axis.title.y = element_text(vjust = -1.5))
gg_assort_gs

# Together
gg_net_mod_betaeps <- grid.arrange(gg_mod_beta, gg_mod_eps, nrow = 1)
ggsave(gg_net_mod, filename = "Output/Networks/NetworkMetrics/Modularity_betaepssweep.svg", 
       height = 23, width = 95, units = "mm")
gg_net_mod <- grid.arrange(gg_mod_beta, gg_mod_eps, gg_mod_gs, nrow = 1)
ggsave(gg_net_mod, filename = "Output/Networks/NetworkMetrics/Modularity_parametersweep.svg", 
       height = 23.1, width = 1.515*95, units = "mm")

gg_net_assort_betaeps <- grid.arrange(gg_assort_beta, gg_assort_eps, nrow = 1)
ggsave(gg_net_assort, filename = "Output/Networks/NetworkMetrics/Assortativity_betaepssweep.svg", 
       height = 23, width = 96, units = "mm")
gg_net_assort <- grid.arrange(gg_assort_beta, gg_assort_eps, gg_assort_gs, nrow = 1)
ggsave(gg_net_assort, filename = "Output/Networks/NetworkMetrics/Assortativity_parametersweep.svg", 
       height = 23.1, width = 1.475*96, units = "mm")
