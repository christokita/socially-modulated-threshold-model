################################################################################
#
# Social interaction model: Analyze parameter space data
#
################################################################################

rm(list = ls())

####################
# Source necessary scripts/libraries
####################
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

####################
# Plot: Beta sweep
####################
load("output/ParameterSpace/GroupSizeBetaSweep_Sigma0-Epsilon0.1.Rdata")
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)

gg_beta <- ggplot(data = entropy, aes(x = n, y = beta, fill = Dind_mean, colour = Dind_mean)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 100, 20), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0.9, 1.25, 0.05), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  # scale_fill_viridis(option = "plasma", direction = -1,
  #                    name = "Behavioral\nspecialization",
  #                    limits = c(0,1)) +
  # scale_colour_viridis(option = "plasma",  direction = -1,
  #                      name = "Behavioral\nspecialization",
  #                      limits = c(0,1)) +
  ylab(expression(beta)) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_beta

if (0.9 %in% unique(entropy$beta)) {
  ggsave(gg_beta, file = "output/ParameterSpace/Plots/BetaGroupSizeSpace_expanded.png", height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_beta, file = "output/ParameterSpace/Plots/BetaGroupSizeSpace_expanded.svg", height = 45, width = 45, units = "mm")
} else {
  ggsave(gg_beta, file = "output/ParameterSpace/Plots/BetaGroupSizeSpace.png", height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_beta, file = "output/ParameterSpace/Plots/BetaGroupSizeSpace.svg", height = 45, width = 45, units = "mm")
}

####################
# Plot: Epsilon sweep
####################
load("output/ParameterSpace/GroupSizeEpsilonSweep_Sigma0-Beta1.1.Rdata")
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)

gg_eps <- ggplot(data = entropy, aes(x = n, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 100, 20), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  ylab(expression(epsilon)) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_eps

if (-0.1 %in% unique(entropy$epsilon)) {
  ggsave(gg_eps, file = "output/ParameterSpace/Plots/EpsilonGroupSizeSpace_expanded.png", height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_eps, file = "output/ParameterSpace/Plots/EpsilonGroupSizeSpace_expanded.svg", height = 45, width = 45, units = "mm")
  
} else {
  ggsave(gg_eps, file = "output/ParameterSpace/Plots/EpsilonGroupSizeSpace.png", height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_eps, file = "output/ParameterSpace/Plots/EpsilonGroupSizeSpace.svg", height = 45, width = 45, units = "mm")
}

