################################################################################
#
# Analyze: The effect of deltas and sigmas on DOL emergence in fixed thresholds
#
################################################################################

rm(list = ls())

source("scripts/util/__Util__MASTER.R")

####################
# Load and bind data
####################
files <- list.files("output/FixedThreshold-SigmaDeltaSweep/Rdata/", pattern = ".Rdata", full.names = T)
for (file in files) {
  load(file)
  if (!exists("entropy")) {
    entropy <- entropies
  } else {
    entropy <- rbind(entropy, entropies)
  }
}

####################
# Summarise
####################
entropy_total <- entropy %>% 
  mutate(sigma = as.factor(sigma)) %>% 
  group_by(sigma, delta, n) %>% 
  summarise(Mean = mean(Dind),
            SE = sd(Dind) / sqrt(length(Dind)), 
            test = length(Dind))

entropy_total$delta <- factor(entropy_total$delta, 
                              levels = c(0.5, 0.6, 0.7, 0.8), 
                              labels = c(expression(paste(delta, "=", 0.5)), 
                                       expression(paste(delta, "=", 0.6)), 
                                       expression(paste(delta, "=", 0.7)), 
                                       expression(paste(delta, "=", 0.8))))

###################
# Plot
####################
library(RColorBrewer)
pal <- brewer.pal(9, "BuPu")[4:9]

gg_deltas <- ggplot(data = entropy_total, aes(x = n, y = Mean, colour = sigma, group = sigma)) +
  geom_line(size = 0.2) +
  # geom_errorbar(aes(ymax = Mean + SE, ymin = Mean - SE)) +
  geom_point(size = 0.6) +
  scale_y_continuous(limits = c(0, 1)) +
  # scale_color_gradient(low = "#9ebcda", high = "#4d004b", limit = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  scale_color_manual(name = expression(paste("Threshold\nvariation (", sigma, ")")),
                     values = pal) +
  theme_classic() +
  xlab("Group size (n)") +
  ylab(expression(paste("Division of labor (", 'D'[indiv], ")"))) +
  facet_grid(~delta, labeller = label_parsed) +
  theme(axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6), 
        legend.key.height = unit(3, "mm"))

gg_deltas
ggsave(gg_deltas, filename = "output/FixedThreshold-SigmaDeltaSweep/FixedThresholdDeltaSigma.png", 
       width = 140, height = 50, units = "mm", dpi = 400)
