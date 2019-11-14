################################################################################
#
# Comparing various specialization plots
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


############### Sweep across beta values ###############
# Compare b=1.1 at long run time to b=1.2 at normal run time. 
# If run time is the main difference in the simulation results, these
# two runs should look similar. n* for b=1.1 is 27.5, for b=1.2

####################
# Load data
####################
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1-LongRun.Rdata")
compiled_data$Model <- "Social_Beta1.1-LongRun"
entropy_data <- compiled_data
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")
compiled_data$Model <- "Social_Beta1.1"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)


load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.2.Rdata")
compiled_data$Model <- "Social_Beta1.2"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)


####################
# Summarise data
####################
# Calculate mean and SE
entropy <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            CI_95 = qnorm(0.975) * sd(Dind) / sqrt(length(n)),
            SD = sd(Dind),
            SE = sd(Dind) / sqrt(length(Dind)))


####################
# Plot
####################
gg_entropy <- ggplot(data = entropy, aes(x = n, colour = Model)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_line(aes(y = Mean),
            size = 0.4) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_color_manual(values = c("#fb9a99", "#e31a1c", "#969696"), 
                     labels = c(expression(paste(italic(beta), " = 1.1")), 
                                expression(paste(italic(beta), " = 1.1, long run")), 
                                expression(paste(italic(beta), " = 1.2"))),
                     name = "Model") +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = c(0.75, 0.2),
        legend.title = element_text(size = 6, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_entropy

ggsave(gg_entropy, file = "output/SpecializationPlots/Sigma0-Epsilon0.1-BetaSweep.svg", 
       height = 45, width = 45, units = "mm")

####################
# Only beta = 1.1
####################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

# Load data
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")
compiled_data$Model <- "Social"
entropy_data <- compiled_data
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0.05-Epsilon0-Beta1.1.Rdata")
compiled_data$Model <- "Fixed"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)

# Calculate mean and SE
entropy <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind),
            SE = sd(Dind) / sqrt(length(Dind)))

# Plot

gg_solo <- ggplot(data = entropy, aes(x = n, colour = Model, fill = Model)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                position = position_dodge(width = 1.5),
                size = 0.3) +
  geom_line(aes(y = Mean),
            position = position_dodge(width = 1.5),
            size = 0.4) +
  geom_point(aes(y = Mean),
             position = position_dodge(width = 1.5),
             size = 1, shape = 21) +
  theme_classic() +
  xlab(expression(paste("Group size (", italic(n), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_color_manual(name = "Threshold",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold",
                    values = c("#ffffff", "#4d4d4d")) +
  scale_linetype_manual(name = "Threshold",
                        values = c("dotted", "solid")) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = "none",
        legend.title = element_text(size = 7, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(5, "mm"),
        axis.ticks = element_line(size = 0.3),
        axis.line = element_line(size = 0.3),
        aspect.ratio = 1)
gg_solo

ggsave(gg_solo, filename = "output/SpecializationPlots/Beta1.1_withfixed.svg", width = 60.5, height = 60.5, units = "mm")
ggsave(gg_solo, filename = "output/SpecializationPlots/Beta1.1withfixed.png", width = 60.5, height = 60.5, units = "mm", dpi = 400)


############### Sweep across epsilon values ###############
rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

####################
# Load data
####################
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.4-Beta1.1.Rdata")
compiled_data$Model <- "Social_Epsilon 0.0"
entropy_data <- compiled_data
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.4-Beta1.1-LongRun.Rdata")
compiled_data$Model <- "Social_Epsilon0.4"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)

####################
# Summarise data
####################
# Calculate mean and SE
entropy <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind),
            SE = sd(Dind) / sqrt(length(Dind)))

####################
# Plot
####################
pal <- brewer.pal(5, "Greens")[c(2, 4, 5)]

gg_entropy <- ggplot(data = entropy, aes(x = n, colour = Model)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_line(aes(y = Mean),
            size = 0.4) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Group size (", italic(n), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_color_manual(values = pal, 
                     labels = c("1.2", "1.1"),
                     name = expression("Interaction bias"(Beta))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = "none",
        legend.title = element_text(size = 7, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(5, "mm"),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_entropy

ggsave(gg_entropy, file = "output/SpecializationPlots/Sigma0-Beta1.1-EpsSweep.png", 
       height = 45, width = 45, units = "mm", dpi = 800)
ggsave(gg_entropy, file = "output/SpecializationPlots/Sigma0-Beta1.1-EpsSweep.svg", 
       height = 48.5, width = 48.5, units = "mm")

