################################################################################
#
# Analyze Networks of Social Interactions
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


filename <- "Sigma001-Eps001-Phi001-ConnectP01-Bias1.1"

# Cutoff for threshold ratio to allow easier plotting
ThreshCutoffValue <- 10

####################
# Compare entropies
####################
# Load social
load("output/Rdata/Sigma001-Eps001-Phi001-ConnectP01-Bias1.1.Rdata")

graphs <- unlist(groups_graphs, recursive = FALSE)
threshMat <- unlist(groups_thresh, recursive = FALSE)
actMat <- unlist(groups_taskDist, recursive = FALSE)

social_graphs <- lapply(1:length(graphs), function(i) {
  # Calculated degree
  degree <- rowSums(graphs[[i]])
  degree <- as.data.frame(degree)
  degree$Id <- row.names(degree)
  # Calculate thresholds
  thresh <- threshMat[[i]] %>% 
    as.data.frame(.) %>% 
    mutate(ThreshRatio = log(Thresh1 / Thresh2),
           Id = row.names(.))
  thresh$ThreshRatio[thresh$ThreshRatio > ThreshCutoffValue] <- Inf
  thresh$ThreshRatio[thresh$ThreshRatio < -ThreshCutoffValue] <- -Inf
  # Calculate actibity
  activity <- actMat[[i]] %>% 
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
load("output/Rdata/Sigma001-FIXED-ConnectP01-Bias1.1.Rdata")

graphs <- unlist(groups_graphs, recursive = FALSE)
threshMat <- unlist(groups_thresh, recursive = FALSE)
actMat <- unlist(groups_taskDist, recursive = FALSE)

fixed_graphs <- lapply(1:length(graphs), function(i) {
  # Calculated degree
  degree <- rowSums(graphs[[i]])
  degree <- as.data.frame(degree)
  degree$Id <- row.names(degree)
  # Calculate thresholds
  thresh <- threshMat[[i]] %>% 
    as.data.frame(.) %>% 
    mutate(ThreshRatio = log(Thresh1 / Thresh2),
           Id = row.names(.))
  thresh$ThreshRatio[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffValue
  thresh$ThreshRatio[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffValue
  # Calculate actibity
  activity <- actMat[[i]] %>% 
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

test <- combined_graphs %>% 
  filter(n == 100)

# Comparing activity and degree, colored by threshold ratio
gg_degree_act <- ggplot(data = social_graphs, 
                        aes(x = ActTotal, y = degree, colour = ThreshRatio, group = n)) +
  geom_point(size = 0.5) +
  theme_bw() +
  scale_colour_gradient2(high = "#d7191c",
                         mid = "#ffffbf", 
                         low = "#2c7bb6", 
                         midpoint = 0, 
                         limits = c(-2, 2),
                         oob = squish) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA)) +
  xlab("Total Activity Level") +
  ylab("Total Degree") +
  facet_wrap( ~ n, scales = "free_y")
gg_degree_act
ggsave(filename = "output/NetworkDataPlots/Social_DegreeVsActivity.png", width = 8, height = 6, units = "in", dpi = 600)


gg_FIX_degree_act <- ggplot(data = fixed_graphs, 
                        aes(x = ActTotal, y = degree, colour = ThreshRatio, group = n)) +
  geom_point() +
  theme_bw() +
  scale_colour_gradient2(high = "#2c7bb6",
                         mid = "#ffffbf", 
                         low = "#d7191c", 
                         midpoint = 0, 
                         limits = c(-2, 2),
                         oob = squish) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA)) +
  facet_wrap( ~ n, scales = "free_y")
gg_FIX_degree_act
ggsave(filename = "output/NetworkDataPlots/Fixed_DegreeVsActivity.png", width = 8, height = 6, units = "in", dpi = 600)


gg_thresh_degree <- ggplot(data = social_graphs, 
                        aes(x = ThreshRatio, y = degree, colour = ActTotal, group = n)) +
  geom_point() +
  theme_bw() +
  scale_colour_gradient2(high = "#6e016b",
                         mid = "#8c6bb1", 
                         low = "#bfd3e6", 
                         limits = c(0, 1),
                         midpoint = 0.5,
                         oob = squish) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA)) +
  facet_wrap( ~ n, scales = "free_y")
gg_thresh_degree


test <- social_graphs %>% filter(n == 100)
qplot(data = test, y = ThreshRatio, x = degree, color = ActTotal) + 
  theme_bw() +
  scale_colour_gradient2(high = "#00441b", mid = "#41ab5d", low = "#c7e9c0", midpoint = 0.5, limits = c(0, 1)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA))

