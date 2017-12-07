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
load("output/Rdata/Sigma0.01-FIXED-Bias1.1.Rdata")

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
                         limits = c(-1, 1),
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
                         limits = c(-1, 1),
                         oob = squish) +
  xlab("Task 1") +
  ylab("Task 2") +
  facet_wrap( ~ n)
gg_activity
ggsave(filename = "output/NetworkDataPlots/Social_TaskPerfVsThreshRatio.png", width = 5, height = 4, units = "in", dpi = 600)

####################
# Compare social network features
####################
# Network "dispersion" (? - standard deviation over mean degree)
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
  summarise(Dispersion = mean(Dispersion))

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
  summarise(Dispersion = mean(Dispersion))

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

# Network homophily
weighted_ratios <- lapply(1:length(soc_graphs), function(i) {
  # Social
  # Get graph and thresh matrix for simulation
  graph <- soc_graphs[[i]]
  thresh <- as.data.frame(soc_threshMat[[i]])
  # Calculate thresh ratio
  threshRatio <- log(thresh$Thresh1 / thresh$Thresh2)
  # Calculate weighted neighbor sum for each individual
  weighted_sum <- graph %*% threshRatio
  # Construct dataframe to return
  to_return <- data.frame(Id = row.names(weighted_sum), 
                          n = length(weighted_sum),
                          ThreshRatio = threshRatio,
                          WeightNeighbor = weighted_sum,
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
  # Construct dataframe to return
  to_return2 <- data.frame(Id = row.names(weighted_sum), 
                          n = length(weighted_sum),
                          ThreshRatio = threshRatio,
                          WeightNeighbor = weighted_sum,
                          Model = "Fixed")
  row.names(to_return) <- NULL
  # Return
  to_return$Model <- as.character(to_return$Model)
  to_return <- rbind(to_return, to_return2)
  return(to_return)
})
weighted_ratios <- do.call("rbind", weighted_ratios)

gg_weighted_ratios <- ggplot(data = weighted_ratios, aes(x = ThreshRatio, y = WeightNeighbor)) +
  geom_point() +
  theme_classic(base_size = 10) +
  theme(aspect.ratio = 1,
        axis.text = element_text(color = "black"))  +
  facet_wrap(~n, scales = "free")
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
