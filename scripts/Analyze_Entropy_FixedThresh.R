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

####################
# Load data
####################
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.01.Rdata")
compiled_data$Model <- "Social_Beta1.01"
entropy_data <- compiled_data
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.05.Rdata")
compiled_data$Model <- "Social_Beta1.05"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")
compiled_data$Model <- "Social_Beta1.1"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.2.Rdata")
compiled_data$Model <- "Social_Beta1.2"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)


load("output/Rdata/_ProcessedData/Entropy/Sigma0.05-Epsilon0-Beta1.1.Rdata")
compiled_data$Model <- "Fixed"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)


####################
# Summarise data
####################
# Calculate mean and SE
entropy <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            SE = sd(Dind) / sqrt(length(Dind)))


####################
# Plot
####################
pal <- c("black", "grey", brewer.pal(5, "Greens")[2:5])

gg_entropy <- ggplot(data = entropy_data, aes(x = n, colour = Model)) +
  geom_point(aes(y = Dind),
             size = 0.8,
             alpha = 0.5,
             stroke = 0) +
  theme_classic() +
  ylab("Division of labor") +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_manual(values = pal, 
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
        aspect.ratio = 1) +
  facet_wrap(~Model)
gg_entropy


gg_entropy <- ggplot(data = entropy, aes(x = n, colour = Model)) +
  geom_line(aes(y = Mean),
            size = 0.4) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  ylab("Division of labor") +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_colour_manual(values = pal, 
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

####################
# Only one data set
####################
look <- entropy %>% 
  filter(Model == "Fixed") %>% 
  mutate(n_log = log10(n))

gg_entropy_check <- ggplot(data = look, aes(x = n_log)) +
  # geom_line(aes(y = Mean),
            # size = 0.4) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  ylab("Division of labor") +
  # scale_x_continuous(breaks = seq(0, 100, 20)) +
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
gg_entropy_check
