################################################################################
#
# Analyze Networks of Social Interactions
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


filename <- "Sigma001-Eps001--Bias1.1"

# Cutoff for threshold ratio to allow easier plotting
ThreshCutoffValue <- 10
ThreshCutoffReplacement <- Inf
ThreshCutoffReplacementColor <- 10

####################
# Compare entropies
####################
# Load social
load("output/Rdata/Sigma0.01-Epsilon0.01-Bias1.1.Rdata")

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
  thresh$ThreshRatio[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffReplacement
  thresh$ThreshRatio[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffReplacement
  thresh$ThreshRatioColor <- thresh$ThreshRatio
  thresh$ThreshRatioColor[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffReplacementColor
  thresh$ThreshRatioColor[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffReplacementColor
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
  thresh$ThreshRatio[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffReplacement
  thresh$ThreshRatio[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffReplacement
  thresh$ThreshRatioColor <- thresh$ThreshRatio
  thresh$ThreshRatioColor[thresh$ThreshRatio > ThreshCutoffValue] <- ThreshCutoffReplacementColor
  thresh$ThreshRatioColor[thresh$ThreshRatio < -ThreshCutoffValue] <- -ThreshCutoffReplacementColor
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

####################
# Compare activity, degree, and thershold ratios
####################
# Comparing activity and degree, colored by threshold ratio
gg_degree_act <- ggplot(data = social_graphs, 
                        aes(x = ActTotal, y = degree, colour = ThreshRatioColor, group = n)) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 10) +
  scale_colour_gradient2(name = "Threshold\nRatio (ln)",
                         high = "#d7191c",
                         mid = "#ffffbf", 
                         low = "#2c7bb6", 
                         midpoint = 0, 
                         limits = c(-1, 1),
                         oob = squish) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        axis.text = element_text(color = "black"),
        strip.background = element_blank()) +
  xlab("Total Activity Level") +
  ylab("Total Degree") +
  scale_x_continuous(breaks = c(0, 0.4, 0.8)) +
  facet_wrap( ~ n, scales = "free_y")
gg_degree_act
ggsave(filename = "output/NetworkDataPlots/Social_DegreeVsActivity.png", width = 5, height = 4, units = "in", dpi = 600)


gg_FIX_degree_act <- ggplot(data = fixed_graphs, 
                        aes(x = ActTotal, y = degree, colour = ThreshRatioColor, group = n)) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 10) +
  scale_colour_gradient2(name = "Threshold\nRatio (ln)",
                         high = "#2c7bb6",
                         mid = "#ffffbf", 
                         low = "#d7191c", 
                         midpoint = 0, 
                         limits = c(-2, 2),
                         oob = squish) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        axis.text = element_text(color = "black"),
        strip.background = element_blank()) +
  facet_wrap( ~ n, scales = "free_y")
gg_FIX_degree_act
ggsave(filename = "output/NetworkDataPlots/Fixed_DegreeVsActivity.png", width = 5, height = 4, units = "in", dpi = 600)

# Comparing threshold ratio and total activity
gg_thresh_degree <- ggplot(data = social_graphs, 
                           aes(x = ThreshRatio, y = degree, colour = ActTotal, group = n)) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 10) +
  scale_colour_gradient2(name = "Total\nActivity",
                         high = "#00441b",
                         mid = "#41ab5d", 
                         low = "#c7e9c0", 
                         midpoint = 0.5, 
                         limits = c(0, 1)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        axis.text = element_text(color = "black"),
        strip.background = element_blank()) +
  xlab("Thresold Ratio") +
  ylab("Degree") +
  facet_wrap( ~ n, scales = "free_y")
gg_thresh_degree
ggsave(filename = "output/NetworkDataPlots/Social_ThreshRatioVsDegree.png", width = 5, height = 4, units = "in", dpi = 600)

gg_FIX_thresh_degree <- ggplot(data = fixed_graphs, 
                           aes(x = ThreshRatio, y = degree, colour = ActTotal, group = n)) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 10) +
  scale_colour_gradient2(name = "Total\nActivity",
                         high = "#00441b",
                         mid = "#41ab5d", 
                         low = "#c7e9c0", 
                         midpoint = 0.5, 
                         limits = c(0, 1)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        axis.text = element_text(color = "black"),
        strip.background = element_blank()) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05)) +
  xlab("Thresold Ratio") +
  ylab("Degree") +
  facet_wrap( ~ n, scales = "free_y")
gg_FIX_thresh_degree
ggsave(filename = "output/NetworkDataPlots/Fixed_ThreshRatioVsDegree.png", width = 5, height = 4, units = "in", dpi = 600)

# Comparing activity and degree, colored by threshold ratio
gg_degree_actRatio <- ggplot(data = social_graphs, 
                        aes(x = ActRatio, y = degree, colour = ThreshRatioColor, group = n)) +
  geom_point(size = 0.5) +
  theme_bw(base_size = 10) +
  scale_colour_gradient2(name = "Threshold\nRatio",
                         high = "#d7191c",
                         mid = "#ffffbf", 
                         low = "#2c7bb6", 
                         midpoint = 0, 
                         limits = c(-1, 1),
                         oob = squish) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA)) +
  xlab("Total Activity Level") +
  ylab("Total Degree") +
  facet_wrap( ~ n, scales = "free_y")
gg_degree_actRatio

# Comparing thresholds
gg_activity <- ggplot(data = social_graphs, 
                             aes(x = Task1, y = Task2, color = ThreshRatioColor, group = n)) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        axis.text = element_text(color = "black"),
        strip.background = element_blank()) +
  scale_colour_gradient2(name = "Threshold\nRatio (ln)",
                         high = "#d7191c",
                         mid = "#ffffbf", 
                         low = "#2c7bb6", 
                         midpoint = 0, 
                         limits = c(-0.5, 0.5),
                         oob = squish) +
  xlab("Task 1") +
  ylab("Task 2") +
  facet_wrap( ~ n)
gg_activity
ggsave(filename = "output/NetworkDataPlots/Social_TaskPerfVsThreshRatio.png", width = 5, height = 4, units = "in", dpi = 600)

# Comparing thresholds

