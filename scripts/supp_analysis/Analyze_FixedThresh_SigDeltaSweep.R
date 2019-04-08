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

entropy_total$delta_label <- factor(entropy_total$delta, 
                              levels = c(0.5, 0.6, 0.7, 0.8), 
                              labels = c(expression(paste(delta, "=", 0.5)), 
                                       expression(paste(delta, "=", 0.6)), 
                                       expression(paste(delta, "=", 0.7)), 
                                       expression(paste(delta, "=", 0.8))))

###################
# Plot
####################
library(RColorBrewer)
library(scales)
pal <- brewer.pal(9, "BuPu")[4:9]


# Line plots
gg_deltas <- ggplot(data = entropy_total, aes(x = n, y = Mean, colour = sigma, group = sigma)) +
  geom_line(size = 0.2) +
  geom_errorbar(aes(ymax = Mean + SE, ymin = Mean - SE), size = 0.2, width = 0) +
  geom_point(size = 0.6) +
  scale_y_continuous(limits = c(0, 1)) +
  # scale_color_gradient(low = "#9ebcda", high = "#4d004b", limit = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  # scale_color_manual(name = expression(paste("Threshold\nvariation (", sigma, ")")),
  #                    values = pal) +
  xlab("Group size (n)") +
  ylab(expression(paste("Division of labor (", 'D'[indiv], ")"))) +
  facet_grid(~delta_label, labeller = label_parsed) +
  theme_ctokita()

gg_deltas
ggsave(gg_deltas, filename = "output/FixedThreshold-SigmaDeltaSweep/FixedThresholdDeltaSigma.png", 
       width = 140, height = 40, units = "mm", dpi = 400)


# heat map
entropy_map <- entropy_total %>% 
  mutate(n = factor(n, levels = unique(entropy_total$n))) %>% 
  filter(delta == 0.8)

pal <- brewer_pal("seq", "Greys")
pal <- pal(9)
pal <- c("#f0f0f0", "#252525")


gg_deltamap <- ggplot(data = entropy_map, aes(x = n, y = sigma, fill = Mean, colour = Mean)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal, name = expression(paste("Division of\nlabor (", 'D'[indiv], ")")),
                         limits = c(0, 0.5)) +
  scale_colour_gradientn(colours = pal,expression(paste("Division of\nlabor (", 'D'[indiv], ")")),
                       limits = c(0, 0.5)) +
  xlab("Group size (n)") +
  scale_x_discrete(breaks = c(5, seq(10, 100, 10))) +
  ylab(expression(paste("Threshold variation (", sigma[j], ")"))) +
  # facet_grid(~delta, labeller = label_parsed) +
  theme_ctokita() +
  theme(aspect.ratio = 1)

gg_deltamap
ggsave(gg_deltamap, filename = "output/FixedThreshold-SigmaDeltaSweep/FixedThresholdDeltaSigma_HeatMap.png", 
       width = 140, height = 40, units = "mm", dpi = 400)
