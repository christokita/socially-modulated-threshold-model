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
load("output/Rdata/_ProcessedData/RankCorr/Sigma0-Epsilon0.1-Beta1.1_RankCorr.Rdata")

####################
# Process data
####################
rank_corr <- compiled_data %>% 
  mutate(Spec = (Task1 + Task2) / 2) %>% 
  group_by(n) %>% 
  summarise(MeanSpec = mean(Spec),
            SESpec = sd(Spec)/sqrt(length(Spec)))


####################
# Plot
####################
gg_spec <- ggplot(data = rank_corr, aes(x = n)) +
  geom_line(aes(y = MeanSpec),
            size = 0.4) +
  geom_errorbar(aes(ymin = MeanSpec - SESpec, ymax = MeanSpec + SESpec),
                width = 0) +
  geom_point(aes(y = MeanSpec),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste("Specialization (rank corr.)"))) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 1)) +
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
gg_spec
