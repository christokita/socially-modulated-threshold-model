################################################################################
#
# Analyze data from long runs of particular parameter combinations
#
################################################################################

##############################  DOL heat map plots  ############################## 
rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

####################
# Load data
####################
# Parameter space data
load("output/ParameterSpace/EpsilonBetaSweep-n80.Rdata")
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)

# Specific parameters (all run for 10x longer, i.e., 500k time steps)
load("output/Rdata/_ProcessedData/_LongSims/Entropy/Sigma0-Epsilon0.5-Beta1.05.Rdata")
long_sims <- entropy_data %>% 
  mutate(epsilon = 0.5, 
         beta = 1.05) %>% 
  group_by(epsilon, beta) %>% 
  summarise(Dind_mean = mean(Dind))

load("output/Rdata/_ProcessedData/_LongSims/Entropy/Sigma0-Epsilon0.55-Beta1.2.Rdata")
entropy_data <- entropy_data %>% 
  mutate(epsilon = 0.55, 
         beta = 1.2) %>% 
  group_by(epsilon, beta) %>% 
  summarise(Dind_mean = mean(Dind))
long_sims <- rbind(long_sims, entropy_data)

load("output/Rdata/_ProcessedData/_LongSims/Entropy/Sigma0-Epsilon0.5-Beta1.1.Rdata")
entropy_data <- entropy_data %>% 
  mutate(epsilon = 0.5, 
         beta = 1.1) %>% 
  group_by(epsilon, beta) %>% 
  summarise(Dind_mean = mean(Dind))
long_sims <- rbind(long_sims, entropy_data)

####################
# Plot
####################
# Graph
gg_betaeps <- ggplot(data = entropy, aes(x = beta, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
  geom_tile() +
  theme_bw() +
  geom_point(data = long_sims, 
             aes(x = beta, y = epsilon, fill = Dind_mean), 
             color = "black",
             shape = 22, 
             size = 2) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  xlab(expression(paste("Interaction Bias ", italic(beta)))) +
  ylab(expression(paste( "Social influence ", italic(epsilon)))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_betaeps

ggsave(gg_betaeps, file = "output/LongSims/BeataEpsSweep_n80_withLongSims.png", height = 45, width = 45, units = "mm", dpi = 400)



##############################  Network property plots  ############################## 

##########################################################
# Modularity
##########################################################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")

p <- 1 #prob of interact
runs <- c("Sigma0-Epsilon0.1_BetaSweep",
          "Sigma0-Beta1.1_EpsSweep")
run_names <- c("Beta", "Epsilon")

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

# Specific parameters (all run for 10x longer, i.e., 500k time steps)
long_runs <- c("Sigma0-Epsilon0.1-Beta1.025",
               "Sigma0-Epsilon0.1-Beta1.05",
               "Sigma0-Epsilon0.1-Beta1.075")

long_run_mod <- lapply(1:length(long_runs), function(run) {
  # Load social networks
  load(paste0("output/Rdata/_ProcessedData/_LongSims/Graphs/", long_runs[run], ".Rdata"))
  # Load threshold matrices
  load(paste0("output/Rdata/_ProcessedData/_LongSims/Thresh/", long_runs[run], ".Rdata"))
  # For each each compute interaction matrix
  # Get graph and make adjacency matrix
  size_graph <- lapply(1:length(graphs_data), function(j) {
    # Format: set diagonal, rescale, and make adj matrix
    this_graph <- graphs_data[[j]]
    diag(this_graph) <- 0
    g <- graph_from_adjacency_matrix(this_graph, mode = "undirected", weighted = TRUE)
    g_clust <- cluster_fast_greedy(g, weights = E(g)$weight)
    # g_membership <- membership(g_clust)
    # mod <- modularity(g, membership = g_membership, weights = E(g)$weight)
    mod <- modularity(g_clust)
    clust_coeff <- transitivity(graph = g, type = "weighted", weights = E(g)$weight)
    # return
    replicate_row <- data.frame(Modularity = mod,
                                ClustCoeff =  mean(clust_coeff, na.rm = TRUE))
    return(replicate_row)
  })
  size_data <- do.call("rbind", size_graph)
  size_data <- size_data %>% 
    mutate(parameter_value = as.numeric(gsub(long_runs[run], pattern = ".*Beta([\\.0-9]+)", replacement = "\\1", perl = T)),
           parameter = "Beta") %>% 
    group_by(parameter_value,
             parameter) %>% 
    summarise(Modul_mean = mean(Modularity),
              Modul_SD = sd(Modularity))
  return(size_data)
})
long_run_mod <- do.call('rbind', long_run_mod)


# Plot
mod_data_beta <- mod_data %>% 
  filter(parameter == "Beta")
gg_mod_beta <- ggplot(mod_data_beta, aes(x = parameter_value, y = Modul_mean, colour = parameter, fill = parameter)) +
  # normal data
  geom_errorbar(aes(ymin = ifelse(Modul_mean - Modul_SD < 0, 0, Modul_mean - Modul_SD), ymax = Modul_mean + Modul_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  # long run data
  geom_errorbar(data = long_run_mod,
                aes(ymin = ifelse(Modul_mean - Modul_SD < 0, 0, Modul_mean - Modul_SD), ymax = Modul_mean + Modul_SD),
                width = 0,
                size = 0.3) +
  geom_point(data = long_run_mod,
             size = 0.8, shape = 21, fill = "red") +
  # Other stuff
  scale_color_manual(name = "Threshold",
                     # values = c("#878787", "#4d4d4d")) +
                     values = c("#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    # values = c("#ffffff", "#4d4d4d")) +
                    values = c("#4d4d4d")) +
  scale_x_continuous(breaks = seq(1, 1.5, 0.05)) +
  scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(-0.0002, 0.031)) +
  xlab(expression(paste("Interaction bias ", italic(beta)))) +
  ylab("Modularity") +
  theme_ctokita() +
  theme(legend.position = "none",
        legend.key.height = unit(0.5, "line"))
gg_mod_beta

ggsave(gg_mod_beta, file = "output/LongSims/modularity_plots.png", width = 45, height = 45, units = "mm", dpi = 400)

#########################################################
# Assortivity
##########################################################
source("scripts/util/__Util__MASTER.R")

p <- 1 #prob of interact
runs <- c("Sigma0-Epsilon0.1_BetaSweep",
          "Sigma0-Beta1.1_EpsSweep")
run_names <- c("Beta", "Epsilon")

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

# Specific parameters (all run for 10x longer, i.e., 500k time steps)
long_runs <- c("Sigma0-Epsilon0.1-Beta1.025",
               "Sigma0-Epsilon0.1-Beta1.05",
               "Sigma0-Epsilon0.1-Beta1.075")

long_run_assort <- lapply(1:length(long_runs), function(run) {
  # Load social networks
  load(paste0("output/Rdata/_ProcessedData/_LongSims/Graphs/", long_runs[run], ".Rdata"))
  # Load threshold matrices
  load(paste0("output/Rdata/_ProcessedData/_LongSims/Thresh/", long_runs[run], ".Rdata"))
  # For each each compute interaction matrix
  # Get graph and make adjacency matrix
  size_graph <- lapply(1:length(graphs_data), function(j) {
    # Get graph and calculate threshold differences
    this_graph <- graphs_data[[j]]
    diag(this_graph) <- 0
    thresh <- as.data.frame(thresh_data[[j]])
    thresh$ThreshBias <- thresh$Thresh1 - thresh$Thresh2
    # Calculate assortmnet
    assort <- assortment.continuous(graph = this_graph, vertex_values = thresh$ThreshBias, weighted = T)
    assort <- assort$r
    to_retun <- data.frame(Assortativity = assort)
    # return
    return(to_retun)
  })
  size_data <- do.call("rbind", size_graph)
  size_data <- size_data %>% 
    mutate(parameter_value = as.numeric(gsub(long_runs[run], pattern = ".*Beta([\\.0-9]+)", replacement = "\\1", perl = T)),
           parameter = "Beta") %>% 
    group_by(parameter_value,
             parameter) %>% 
    summarise(Assort_mean = mean(Assortativity),
              Assort_SD = sd(Assortativity))
  return(size_data)
})
long_run_assort <- do.call('rbind', long_run_assort)

# Plot
assort_data_beta <- assort_data %>% 
  filter(parameter == "Beta")
gg_assort_beta <- ggplot(data = assort_data_beta, aes(x = parameter_value, y = Assort_mean,
                                            colour = parameter, group = parameter, fill = parameter)) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dotted") +
  # Regular data
  geom_errorbar(aes(ymin = Assort_mean - Assort_SD, ymax = Assort_mean + Assort_SD),
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8, shape = 21) +
  # long run data
  geom_errorbar(data = long_run_assort,
                aes(ymin = Assort_mean - Assort_SD, ymax = Assort_mean + Assort_SD),
                width = 0,
                size = 0.3) +
  geom_point(data = long_run_assort,
             size = 0.8, shape = 21, fill = "red") +
  # Other stuff
  scale_color_manual(name = "Threshold",
                     # values = c("#878787", "#4d4d4d")) +
                     values = c("#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    # values = c("#ffffff", "#4d4d4d")) +
                    values = c("#4d4d4d")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(-0.04, 0.1, 0.02), limits = c(-0.02, 0.062)) +
  xlab(expression(paste("Interaction bias ", italic(beta)))) +
  ylab("Assortativity") +
  theme_ctokita() +
  theme(legend.position = "none",
        axis.title.y = element_text(vjust = -1.5))
gg_assort_beta

ggsave(gg_assort_beta, file = "output/LongSims/assortativity_plots.png", width = 45, height = 45, units = "mm", dpi = 400)

