################################################################################
#
# Social interaction model: check results of simulation length test
#
################################################################################

rm(list = ls())

####################
# Source necessary scripts/libraries
####################
source("scripts/util/__Util__MASTER.R")
library(scales)

####################
# Load data
####################
path <- "output/SimLength/CheckSimLength_Thresh50_Sigma0-Epsilon0.05-Beta1.1.Rdata"
load(path)

file_name <- gsub("^.*(Thresh[\\. 0-9]*.*)\\.Rdata$", "\\1", path, perl = T)

########## Timeseries data #########

####################
# Summarise
####################
sim_data <- as.data.frame(parallel_simulations)
sim_data <- sim_data %>% 
  select(-sim, -chunk) %>% 
  group_by(n, t) %>% 
  summarise_all(funs(mean))
sim_data$n <- factor(sim_data$n, levels = c(5, 10, 20, 30, 40, 50, 100))
sim_data <- sim_data %>% 
  mutate(TotalStimAvg = Stim1Avg + Stim2Avg,
         Thresh1Gap = Thresh1Max - Thresh1Min,
         Thresh2Gap = Thresh2Max - Thresh2Min) %>% 
  mutate(ThreshGap = (Thresh1Gap + Thresh2Gap) / 2,
         ThreshMax = (Thresh1Max + Thresh2Max) / 2,
         ThreshMin = (Thresh1Min + Thresh2Min) / 2,
         ThreshSD = (Thresh1SD + Thresh2SD) / 2)


####################
# Plot - DOL
####################
pal <- brewer_pal(11, palette = "YlGnBu")
cols <- pal(n = length(unique(sim_data$n)) + 1)[2:(length(unique(sim_data$n))+1)]

ggplot(data = sim_data, aes(x = t, y = Dind, group = n, col = n)) +
  geom_vline(xintercept = 20000, linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = 40000, linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = 50000, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 60000, linetype = "dotted", alpha = 0.3) +
  geom_line() +
  theme_classic() +
  scale_color_manual(values = cols) +
  ylab("Division of Labor") +
  theme(aspect.ratio = 1)

ggsave(file = paste0("output/SimLength/DOL/", file_name, ".png"), scale = 0.6, dpi = 300)

####################
# Plot - Stimulus Levels
####################
pal <- brewer_pal(11, palette = "YlGnBu")
cols <- pal(n = length(unique(sim_data$n)) + 1)[2:(length(unique(sim_data$n))+1)]

ggplot(data = sim_data, aes(x = t, y = TotalStimAvg, group = n, col = n)) +
  geom_vline(xintercept = 20000, linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = 40000, linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = 50000, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 60000, linetype = "dotted", alpha = 0.3) +
  geom_line() +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(aspect.ratio = 1) +
  ylab("Average Stimulus Level")

ggsave(file = paste0("output/SimLength/StimulusLevel/", file_name, ".png"), scale = 0.6, dpi = 300)


####################
# Plot - Thresh SD
####################
pal <- brewer_pal(11, palette = "YlGnBu")
cols <- pal(n = length(unique(sim_data$n)) + 1)[2:(length(unique(sim_data$n))+1)]

ggplot(data = sim_data, aes(x = t, y = ThreshSD, group = n, col = n)) +
  geom_vline(xintercept = 20000, linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = 40000, linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = 50000, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 60000, linetype = "dotted", alpha = 0.3) +
  geom_line() +
  theme_classic() +
  ylab("Threshold Variation") +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = cols)

ggsave(file = paste0("output/SimLength/ThresholdVariation/", file_name, ".png"), scale = 0.6, dpi = 300)

####################
# Plot - Thresh Gap
####################
pal <- brewer_pal(11, palette = "YlGnBu")
cols <- pal(n = length(unique(sim_data$n)) + 1)[2:(length(unique(sim_data$n))+1)]

ggplot(data = sim_data, aes(x = t, y = ThreshGap, group = n, col = n)) +
  geom_vline(xintercept = 20000, linetype = "dotted") +
  geom_vline(xintercept = 50000, linetype = "dashed") +
  geom_vline(xintercept = 20000, linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = 40000, linetype = "dotted", alpha = 0.3) +
  geom_vline(xintercept = 50000, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 60000, linetype = "dotted", alpha = 0.3) +
  geom_line() +
  theme_classic() +
  ylab("Threshold Gap") +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = cols)

ggsave(file = paste0("output/SimLength/ThresholdGap/", file_name, ".png"), scale = 0.6, dpi = 300)


########## Distribution #########

####################
# Plot distribution
####################
sim_data <- as.data.frame(parallel_simulations)
sim_data <- sim_data %>% 
  select(-sim, -chunk)
sim_data$n <- factor(sim_data$n, levels = c(5, 10, 20, 30, 40, 50, 100))
sim_data <- sim_data %>% 
  mutate(TotalStimAvg = Stim1Avg + Stim2Avg,
         Thresh1Gap = Thresh1Max - Thresh1Min,
         Thresh2Gap = Thresh2Max - Thresh2Min) %>% 
  mutate(ThreshGap = (Thresh1Gap + Thresh2Gap) / 2,
         ThreshMax = (Thresh1Max + Thresh2Max) / 2,
         ThreshMin = (Thresh1Min + Thresh2Min) / 2,
         ThreshSD = (Thresh1SD + Thresh2SD) / 2) %>% 
  filter(t %in% c(20000, 30000, 40000, 50000, 60000, 70000))

pal <- brewer_pal(11, palette = "YlGnBu")
cols <- pal(n = length(unique(sim_data$n)) + 1)[2:(length(unique(sim_data$n))+1)]

ggplot(data = sim_data, aes(x = t, y = ThreshGap, group = n, col = n)) +
  geom_point(position = position_jitterdodge()) +
  theme_classic() +
  scale_color_manual(values = cols)


ggplot(data = sim_data, aes(x = ThreshGap, group = n, fill = n)) +
  geom_histogram() +
  theme_classic() +
  scale_fill_manual(values = cols) +
  facet_grid(t ~ n)






                