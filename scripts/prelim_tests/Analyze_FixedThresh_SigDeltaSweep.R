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
files <- list.files("output/FixedThreshold-SigmaDeltaSweep/Rdata/", full.names = T)
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
gg_deltas <- ggplot(data = entropy_total, aes(x = n, y = Mean, colour = sigma, group = sigma)) +
  geom_line(position = position_dodge(width = 5)) +
  geom_errorbar(aes(ymax = Mean + SE, ymin = Mean - SE),
                position = position_dodge(width = 5)) +
  geom_point(position = position_dodge(width = 5)) +
  scale_y_continuous(limits = c(0, 1)) +
  # scale_color_gradient(low = "#9ebcda", high = "#4d004b", limit = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  scale_color_brewer(name = expression(paste("Threshold\nvariation (", sigma, ")")),
                     palette = "RdYlBu",
                     direction = -1) +
  theme_classic() +
  xlab("Group size (n)") +
  ylab(expression(paste("Division of labor (", 'D'[indiv], ")"))) +
  facet_grid(~delta, labeller = label_parsed)

gg_deltas
ggsave(gg_deltas, filename = "output/FixedThreshold-SigmaDeltaSweep/plots")
