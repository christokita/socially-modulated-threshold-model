################################################################################
#
# DOL group size plots: social thresholds vs. fixed thresholds vs. social threshold with initial variation
#
################################################################################


rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

####################
# Load data
####################
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")
compiled_data$Model <- "Social"
entropy_data <- compiled_data
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0.05-Epsilon0.1-Beta1.1.Rdata")
compiled_data$Model <- "Social with variation"
entropy_data <- rbind(entropy_data, compiled_data)
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

gg_comparison <- ggplot(data = entropy, aes(x = n, colour = Model, fill = Model)) +
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
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20))) +
  scale_color_manual(name = "Threshold",
                     values = c("#878787", "#4d4d4d", "#ef8a62")) +
  scale_fill_manual(name = "Threshold",
                    values = c("#ffffff", "#4d4d4d", "#ef8a62")) +
  scale_linetype_manual(name = "Threshold",
                        values = c("dotted", "solid", "solid")) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = c(0.7, 0.7),
        legend.title = element_text(size = 7, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(5, "mm"),
        axis.ticks = element_line(size = 0.3),
        axis.line = element_line(size = 0.3),
        aspect.ratio = 1)
gg_comparison

ggsave(gg_comparison, filename = "output/SpecializationPlots/SocialThresholds_withvariation.png", width = 60.5, height = 60.5, units = "mm", dpi = 400)

