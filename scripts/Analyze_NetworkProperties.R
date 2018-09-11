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
# Percentage of non-random interactions
####################
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
    #  Calcualte 99% CI interval of interaction rate
    edgelist_sig <- all_edgelist %>% 
      group_by(From, To, Interaction) %>% 
      # filter(!is.na(Weight)) %>% 
      summarise(samp_mean = mean(Weight),
                samp_sd = sd(Weight),
                samples = length(Weight)) %>% 
      mutate(error = qt(0.995,df = samples-1) * samp_sd/sqrt(samples),
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
  geom_line(size = 0.4, aes(linetype = Model)) +
  geom_point(size = 0.8, shape = 21) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(name = "Threshold type",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold type",
                    values = c("#ffffff", "#4d4d4d")) +
  scale_linetype_manual(name = "Threshold type",
                        values = c("dotted", "solid")) +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab("% Non-random interactions") +
  theme_ctokita() +
  theme(aspect.ratio = 1,
        legend.position = "none")

gg_interactions
ggsave(gg_interactions, filename = "Output/Networks/NetworkMetrics/PercentNonRandomInteractions.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_interactions, filename = "Output/Networks/NetworkMetrics/PercentNonRandomInteractions.svg", 
       height = 45, width = 45, units = "mm")


####################
# Assortivity
####################
network_correlations <- lapply(1:length(runs), function(run) {
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
      # Get graph and calculate threshold differences
      this_graph <- graphs[[j]]
      diag(this_graph) <- NA
      thresh <- as.data.frame(thresh_data[[i]][j])
      thresh$ThreshDiff <- thresh$Thresh2 - thresh$Thresh1
      thresh$ThreshRatio <- log(thresh$Thresh1 / thresh$Thresh2)
      # Multiply
      social_interaction <- thresh$ThreshDiff * this_graph
      social_interaction <- t(social_interaction)
      effective_interactions <- rowSums(social_interaction, na.rm = T)
      to_retun <- data.frame(n = nrow(this_graph), Correlation = cor(effective_interactions, thresh$ThreshDiff))
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
correlation_data <- do.call('rbind', network_correlations)
correlation_data <- correlation_data %>% 
  group_by(Model, n) %>% 
  summarise(Corr_mean = mean(Correlation),
            Corr_SE = sd(Correlation)/length(Correlation))

# Plot
gg_correlation <- ggplot(data = correlation_data, aes(x = n, y = Corr_mean, 
                                                      colour = Model, group = Model, fill = Model)) +
  geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = Corr_mean - Corr_SE, ymax = Corr_mean + Corr_SE),
                width = 0) +
  geom_point(size = 0.8, shape = 21) +
  scale_color_manual(name = "Threshold type",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold type",
                    values = c("#ffffff", "#4d4d4d")) +
  scale_linetype_manual(name = "Threshold type",
                        values = c("dotted", "solid")) +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab("Social netowrk correlation") +
  theme_ctokita() +
  theme(aspect.ratio = 1,
        legend.position = c(0.8, 0.8))

gg_correlation
ggsave(gg_correlation, filename = "Output/Networks/NetworkMetrics/CorrelationInNetwork.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_correlation, filename = "Output/Networks/NetworkMetrics/CorrelationInNetwork.svg", 
       height = 45, width = 45, units = "mm")
