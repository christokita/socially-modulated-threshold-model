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

filename <- "Sigma0.05-Eps0.01--Bias1.1"

# Cutoff for threshold ratio to allow easier plotting
ThreshCutoffValue <- 10
ThreshCutoffReplacement <- Inf
ThreshCutoffReplacementColor <- 10

####################
# Load data
####################
# Load social
load("output/Rdata/Sigma0.05-Epsilon0.01-Bias1.1.Rdata")

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
           Id = row.names(.))
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
           Id = row.names(.))
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

test <- weighted_ratios %>% filter(Model == "Fixed")

gg_weighted_ratios <- ggplot(data = test, aes(x = ThreshRatio, y = WeightNeighbor)) +
  geom_point() +
  theme_classic(base_size = 10) +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"))  +
  facet_wrap(~n)
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
gg_disparity <- ggplot(data = weighted_ratios, aes(x = ThreshRatio, y = Disparity)) +
  geom_point(size = 0.2) +
  theme_classic(base_size = 10) +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"))  +
  facet_wrap(~n, scales = "free")
gg_disparity


###### Modularity ###### 
