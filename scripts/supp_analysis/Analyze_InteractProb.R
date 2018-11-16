################################################################################
#
# Analyzing the effect of the interaction parameter
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)
library(viridis)
library(ggridges)


####################
# Load and process data
####################
# Load entropy
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1-P0.5.Rdata")
entropy_data <- compiled_data 
entropy_data$Model <- "Ep0.1-P0.5"

load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")
compiled_data$Model <- "Ep0.1-P1"
entropy_data <- rbind(entropy_data, compiled_data)

# Summarise
entropy_data <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            SE = sd(Dind) / sqrt(length(Dind)))

####################
# Plot
####################
pal <- c("#a6cee3", "#1f78b4", "#33a02c")

gg_int <- ggplot(data = entropy_data, aes(x = n, colour = Model)) +
  geom_line(aes(y = Mean),
            size = 0.4) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_color_manual(values = pal, 
                     labels = c(expression(paste(epsilon, " = 0.1, ", p, " = 0.5")),
                                expression(paste(epsilon, " = 0.1, ", p, " = 1"))),
                     name = "Model") +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = "right",
        legend.title = element_text(size = 7, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(5, "mm"),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        legend.text.align = 0,
        aspect.ratio = 1)
gg_int

ggsave("output/InteractionProbability/EntropyComparison.png", width = 80, heigh = 45, units = "mm", dpi = 400)


