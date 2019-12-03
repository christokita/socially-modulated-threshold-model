################################################################################
#
# Testing whether the threshold limit is causing artifacts for high social influence
#
################################################################################
#
# Here we are testing wehther the decrease in DOL at high levels of social influence
# is due to the upper limit on threhsold values. In this simulation, we put the upper bound
# for threshold at 1000, instead of the normal 100.

rm(list = ls())
source("scripts/util/__Util__MASTER.R")


# ------------------------------ DOL by epsilon - LONG RUN----------------------
####################
# Load and process data
####################
load('output/Rdata/_ProcessedData/Entropy/Sigma0-Beta1.1_EpsSweep-NoThreshLimit-LongRun.Rdata')
high_thresh <- compiled_data %>% 
  mutate(Model = "high_thresh") %>% 
  group_by(Model, epsilon, beta) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind)) %>% 
  as.data.frame() %>% 
  mutate(epsilon = round(epsilon, digits = 3))


load('output/Rdata/_ProcessedData/Entropy/Sigma0-Beta1.1_EpsSweep-LongRun.Rdata')
normal_thresh <- compiled_data %>% 
  mutate(Model = "normal_thresh") %>% 
  group_by(Model, epsilon, beta) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind)) %>% 
  as.data.frame() %>% 
  mutate(epsilon = round(epsilon, digits = 3))

entropy_data <- rbind(high_thresh, normal_thresh)


####################
# Plot entropy plots
####################
gg_comp_long <- ggplot(entropy_data, aes(x = epsilon, y = Mean, group = Model, color = Model)) +
  # geom_vline(xintercept = 0.5659957, size = 0.3, linetype = "dashed") + #analytical result for epsilon* 
  geom_vline(xintercept = 0.4969125, size = 0.3, linetype = "dashed", color = "grey40") + #analytical result for epsilon*
  geom_vline(xintercept = 0, size = 0.3, linetype = "dashed", color = "grey40") + #analytical result for epsilon*
  geom_errorbar(aes(ymin = ifelse((Mean - SD) > 0, Mean - SD, 0), ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Social influence ", italic(epsilon)))) +
  ylab(expression(paste("Division of labor ", italic(D[indiv])))) +
  scale_color_manual(name = "Thresh. limits",
                     values = c("#a6cee3", "#1f78b4"),
                     labels = c(expression(paste("[0, ", infinity, ")")), "[0, 100]")) +
  # ggtitle(expression(paste(italic(epsilon), "= 0.4, ", italic(beta), "= 1.1"))) +
  scale_x_continuous(breaks = seq(0, 0.6, 0.1)) +
  theme(title = element_text(size = 6),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = c(0.3, 0.2),
        legend.title = element_text(size = 6, 
                                    face = "bold", 
                                    hjust = 5),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        legend.background = element_blank(),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_comp_long

ggsave(gg_comp_long, file = "output/SpecializationPlots/ThresholdLimitComparison_long.png", height = 45, width = 45, units = "mm")
ggsave(gg_comp_long, file = "output/SpecializationPlots/ThresholdLimitComparison_long.svg", height = 45, width = 45, units = "mm")


# ------------------------------ DOL by epsilon----------------------
####################
# Load and process data
####################
load('output/Rdata/_ProcessedData/Entropy/Sigma0-Beta1.1_EpsSweep-NoThreshLimit.Rdata')
high_thresh <- compiled_data %>% 
  mutate(Model = "high_thresh") %>% 
  group_by(Model, epsilon, beta) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind)) %>% 
  as.data.frame() %>% 
  mutate(epsilon = round(epsilon, digits = 3))


load("output/ParameterSpace/EpsilonBetaSweep-n80.Rdata")
normal_thresh <- entropy %>% 
  mutate(Model = "normal_thresh") %>% 
  filter(beta == 1.1) %>% 
  select(Model, epsilon, beta, Dind_mean, Dind_SD)
names(normal_thresh) <- c("Model", "epsilon", "beta", "Mean", "SD")

eps_values <- unique(normal_thresh$epsilon)

entropy_data <- rbind(high_thresh, normal_thresh)
entropy_data <- entropy_data %>%
  filter(epsilon %in% eps_values)

####################
# Plot entropy plots
####################
gg_comp <- ggplot(entropy_data, aes(x = epsilon, y = Mean, group = Model, color = Model)) +
  # geom_vline(xintercept = 0.5659957, size = 0.3, linetype = "dashed") + #analytical result for epsilon* 
  geom_vline(xintercept = 0.4969125, size = 0.3, linetype = "dashed", color = "grey40") + #analytical result for epsilon*
  geom_vline(xintercept = 0, size = 0.3, linetype = "dashed", color = "grey40") + #analytical result for epsilon*
  geom_errorbar(aes(ymin = ifelse((Mean - SD) > 0, Mean - SD, 0), ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Social influence ", italic(epsilon)))) +
  ylab(expression(paste("Division of labor ", italic(D[indiv])))) +
  scale_color_manual(name = "Thresh. limits",
                     values = c("#a6cee3", "#1f78b4"),
                     labels = c(expression(paste("[0, ", infinity, ")")), "[0, 100]")) +
  # ggtitle(expression(paste(italic(epsilon), "= 0.4, ", italic(beta), "= 1.1"))) +
  scale_x_continuous(breaks = seq(0, 0.6, 0.1)) +
  theme(title = element_text(size = 6),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = c(0.3, 0.2),
        legend.title = element_text(size = 6, 
                                    face = "bold", 
                                    hjust = 5),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        legend.background = element_blank(),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_comp

ggsave(gg_comp, file = "output/SpecializationPlots/ThresholdLimitComparison.png", height = 45, width = 45, units = "mm")
ggsave(gg_comp, file = "output/SpecializationPlots/ThresholdLimitComparison.svg", height = 45, width = 45, units = "mm")


# ------------------------------ DOL by beta----------------------
####################
# Load and process data
####################
load('output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.4_BetaSweep-HighThreshLimit.Rdata')
high_thresh <- compiled_data %>% 
  mutate(Model = "high_thresh") %>% 
  group_by(Model, epsilon, beta) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind)) %>% 
  as.data.frame()


load("output/ParameterSpace/EpsilonBetaSweep-n80.Rdata")
normal_thresh <- entropy %>% 
  mutate(Model = "normal_thresh") %>% 
  filter(epsilon == 0.4) %>% 
  select(Model, epsilon, beta, Dind_mean, Dind_SD)
names(normal_thresh) <- c("Model", "epsilon", "beta", "Mean", "SD")

entropy_data <- rbind(high_thresh, normal_thresh)

####################
# Plot entropy plots
####################
gg_comp <- ggplot(entropy_data, aes(x = beta, y = Mean, group = Model, color = Model)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Group Size ", italic(n)))) +
  ylab(expression(paste("Division of labor ", italic(D[indiv])))) +
  scale_color_manual(name = "Thresh. limits",
                     values = c("#a6cee3", "#1f78b4"),
                     labels = c("[0, 1,000]", "[0, 100]")) +
  # ggtitle(expression(paste(italic(epsilon), "= 0.4, ", italic(beta), "= 1.1"))) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  theme(title = element_text(size = 6),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = c(0.8, 0.2),
        legend.title = element_text(size = 6, 
                                    face = "bold", 
                                    hjust = 5),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        legend.background = element_blank(),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_comp

ggsave(gg_comp, file = "output/SpecializationPlots/ThresholdLimitComparison.png", height = 45, width = 45, units = "mm")
ggsave(gg_comp, file = "output/SpecializationPlots/ThresholdLimitComparison.svg", height = 45, width = 45, units = "mm")


# ------------------------------ DOL by group size ------------------------------
####################
# Load and process data
####################
load('output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.4-Beta1.1-HighThreshLimit.Rdata')
high_thresh <- compiled_data
high_thresh$Model <- "high_thresh"

load('output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.4-Beta1.1.Rdata')
normal_thresh <- compiled_data
normal_thresh$Model <- "normal_thresh"

entropy_data <- rbind(high_thresh, normal_thresh)
entropy_data <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

####################
# Plot entropy plots
####################
gg_comp <- ggplot(entropy_data, aes(x = n, y = Mean, group = Model, color = Model)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_line(aes(y = Mean),
            size = 0.4) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Group Size ", italic(n)))) +
  ylab(expression(paste("Division of labor ", italic(D[indiv])))) +
  scale_color_manual(name = "Thresh. limits",
                     values = c("#a6cee3", "#1f78b4"),
                     labels = c("[0, 1,000]", "[0, 100]")) +
  # ggtitle(expression(paste(italic(epsilon), "= 0.4, ", italic(beta), "= 1.1"))) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme(title = element_text(size = 6),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = c(0.8, 0.2),
        legend.title = element_text(size = 6, 
                                    face = "bold", 
                                    hjust = 5),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        legend.background = element_blank(),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_comp
ggsave(gg_comp, file = "output/SpecializationPlots/ThresholdLimitComparison.png", height = 45, width = 45, units = "mm")
ggsave(gg_comp, file = "output/SpecializationPlots/ThresholdLimitComparison.svg", height = 45, width = 45, units = "mm")
  

# ------------------------------ Thresholds ------------------------------
####################
# Load and process data
####################
files <- list.files('output/Rdata/_ProcessedData/Thresh/Sigma0-Epsilon0.4-Beta1.1-HighThreshLimit/', full.names = T)
group_thresh <- lapply(files, function(file) {
  
  # Load
  load(file)
  
  # Bind
  thresh_data <- as.data.frame(do.call("rbind", listed_data))
  
  # # Summarise
  # thresh_sum <- thresh_data %>% 
  #   mutate(run = paste0(sim, "-", chunk)) %>% 
  #   group_by(n, run) %>% 
  #   summarise(Thresh1SD = sd(Thresh1),
  #             Thresh2SD = sd(Thresh2),
  #             ThreshSD = (sd(Thresh1)  + sd(Thresh2)) / 2) %>% 
  #   select(n, ThreshSD)
  return(thresh_data)
})
# Bind 
all_thresh <- do.call('rbind', group_thresh)


####################
# Plot distirubtion
####################
library(ggridges)
library(viridis)
# Plot
gg_threshvar <- ggplot(data = all_thresh, 
                       aes(x = Thresh1, y = n, fill = n, group = n)) +
  theme_invisible() +
  geom_density_ridges(size = 0.1, stat = "binline", bins = 100) +
  xlab(expression(paste("Task 1 threshold ", italic(theta[i1])))) +
  ylab(expression(paste("Group Size ", italic(n)))) +
  # scale_x_continuous(breaks = seq(0, 100, 25), 
  #                    # limits = c(0, 1),
  #                    expand = c(0.03, 0)) +
  scale_y_continuous(breaks = c(5, seq(10, 50, 10)),
                     expand = c(0.03, 0)) +
  scale_fill_viridis() +
  scale_color_viridis() +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, color = "black"),
        axis.title.y = element_text(size = 7, color = "black", face = 'italic'),
        legend.position = "none",
        strip.text = element_text(size = 6, face = "bold"),
        axis.ticks = element_line(size = 0.2, color = "black"))

gg_threshvar


# ------------------------------ Thresholds over time  ------------------------------
####################
# Load normal threshold limit (epsilon = 0.4, beta = 1.1)
####################
library(scales)
# Normal threshold limit
load("output/ThresholdTime/ThresholdLimits/ThreshMax-100_Eps0.4.Rdata")

gg_threshtime_100 <- ggplot(thresh_time, aes(x = t, y = Threshold, group = Id)) +
  geom_line(size = 0.1, alpha = 0.1, colour = "#1f78b4") +
  scale_x_continuous(name = expression(paste("Time step ", italic(t))),
                     breaks = seq(0, 50000, 10000),
                     labels = c("", "10,000", "", "30,000", "", "50,000"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = expression(paste("Task 1 threshold ", italic(theta[i1,t]))),
                     limits = c(0, 100),
                     breaks = seq(0, 100, 50)) +
  theme_ctokita() +
  theme(#axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 0.7))

gg_threshtime_100
ggsave(gg_threshtime_100, file = "output/ThresholdTime/ThresholdLimits/Max-100.png", height = 26, width = 40, units = "mm", dpi = 500)
ggsave(gg_threshtime_100, file = "output/ThresholdTime/ThresholdLimits/ExampleforPPT.png", height = 34, width = 45, units = "mm", dpi = 500)


# High threshold limit
load("output/ThresholdTime/ThresholdLimits/ThreshMax-Inf_Eps0.4.Rdata")

gg_threshtime_1000 <- ggplot(thresh_time, aes(x = t, y = Threshold, group = Id)) +
  geom_line(size = 0.1, alpha = 0.15, colour = "#8bbeda") +
  scale_x_continuous(name = expression(paste("Time step ", italic(t))),
                     breaks = seq(0, 50000, 10000),
                     labels = c("", "10,000", "", "30,000", "", "50,000"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = expression(paste("Task 1 threshold ", italic(theta[i1,t]))),
                     limits = c(0, 1200),
                     breaks = seq(0, 1200, 600),
                     label = comma) +
  theme_ctokita() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 0.7))

gg_threshtime_1000
ggsave(gg_threshtime_1000, file = "output/ThresholdTime/ThresholdLimits/Max-Inf.png", height = 26, width = 42, units = "mm", dpi = 500)
