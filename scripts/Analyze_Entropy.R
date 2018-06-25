################################################################################
#
# Comparing various specialization plots
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
compiled_data$Model <- "Social_Beta1.1"
entropy_data <- compiled_data
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.2.Rdata")
compiled_data$Model <- "Social_Beta"
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
gg_entropy <- ggplot(data = entropy, aes(x = n, colour = Model)) +
  geom_line(aes(y = Mean),
            size = 0.4) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0) +
  geom_point(aes(y = Mean),
             size = 1) +
  theme_classic() +
  ylab("Division of Labor") +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_color_manual(values = c("#41ab5d", "#005a32"), 
                     labels = c("1.2", "1.1"),
                     name = expression("Interaction bias"(Beta))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.position = "right",
        legend.title = element_text(size = 7, 
                                    face = "bold"),
        legend.text = element_text(size = 6,),
        # legend.key.height = unit(5, "mm"),
        legend.key.width = unit(5, "mm"),
        axis.ticks = element_line(size = 0.3),
        axis.line = element_line(size = 0.3),
        aspect.ratio = 1)
gg_entropy

ggsave(gg_entropy, file = "output/SpecializationPlots/Sigma0-Epsilon0.1-Beta1.1.png", 
       height = 45, width = 45, units = "mm", dpi = 800)
