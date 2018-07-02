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
# Load data
####################
load("output/ParameterSpace/GroupSizeEpsilonSweep_Sigma0-Beta1.1.Rdata")

####################
# Plot: Beta sweep
####################
pal <- brewer_pal("seq", "YlGnBu")
pal <- pal(9)

gg_beta <- ggplot(data = entropy, aes(x = n, y = beta, fill = Dind_mean, colour = Dind_mean)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 100, 20), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(1, 1.25, 0.05), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  ylab(expression(beta)) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        axis.ticks = element_line(size = 0.3),
        panel.background = element_rect(fill = NA, size = 0.3),
        aspect.ratio = 1)
gg_beta

ggsave(gg_beta, file = "output/ParameterSpace/Plots/BetaGroupSizeSpace.png", height = 45, units = "mm", dpi = 400)


####################
# Plot: Epsilon sweep
####################
pal <- brewer_pal("seq", "YlGnBu")
pal <- pal(9)

gg_eps <- ggplot(data = entropy, aes(x = n, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 100, 20), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), 
                     expand = c(0,0),
                     limits = c(0, 0.6)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  ylab(expression(epsilon)) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        axis.ticks = element_line(size = 0.3),
        panel.background = element_rect(fill = NA, size = 0.3),
        aspect.ratio = 1)
gg_eps

ggsave(gg_eps, file = "output/ParameterSpace/Plots/EpsilonGroupSizeSpace.png", height = 45, units = "mm", dpi = 400)

