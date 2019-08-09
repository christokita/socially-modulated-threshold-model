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


# Beta sweeps -----------------------------------------------------------------

####################
# Plot: Beta sweep - positive social influence
####################
load("output/ParameterSpace/GroupSizeBetaSweep_Sigma0-Epsilon0.1.Rdata")
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)

# library(viridis)
# pal <- viridis(9, option = "plasma")

# As separate plots
for (bias in c("Homophily", "Heterophily")) {
  # Filter
  if (bias == "Homophily") {
    entropy_filt <- entropy %>% 
      filter(beta >= 1)
  } else {
    entropy_filt <- entropy %>% 
      filter(beta <= 1)
  }
  # Graph
  gg_beta <- ggplot(data = entropy_filt, aes(x = n, y = beta, fill = Dind_mean, colour = Dind_mean)) +
    geom_tile() +
    theme_bw() +
    scale_x_continuous(breaks = c(5, seq(20, 100, 20)), 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0.75, 1.25, 0.05), 
                       expand = c(0,0)) +
    scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
    scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                           limits = c(0, 1)) +
    xlab(expression(paste("Group Size (", italic(n), ")"))) +
    ylab(expression(paste( "Interaction bias (", italic(beta), ")"))) +
    theme(axis.text = element_text(colour = "black", size = 6),
          axis.title = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.height = unit(5, "mm"),
          legend.key.width = unit(2, "mm"),
          legend.position = "none",
          axis.ticks = element_line(size = 0.3, color = "black"),
          panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
          aspect.ratio = 1)
  gg_beta
  # Save
  file_png <- paste0("output/ParameterSpace/Plots/BetaSweep_Positive", bias,  ".png")
  file_svg <- paste0("output/ParameterSpace/Plots/svg/BetaSweep_Positive", bias,  ".svg")
  ggsave(gg_beta, file = file_png, height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_beta, file = file_svg, height = 45, width = 45, units = "mm")
}

# Analytical results for group size above which full DOL should emerge (n*)
# Not accounting for double interact
analytical_data <- data.frame(beta = seq(1.001, 1.255, 0.0001), n = rep(NA, length(seq(1.001, 1.255, 0.0001))))
for (i in 1:nrow(analytical_data)) {
  analytical_data[i, 2] <- (2 * analytical_data[i, 1]) / (0.8 * (analytical_data[i, 1] - 1))
}
analytical_data <- analytical_data %>% 
  filter(n < 102.5)

# As one plot
gg_beta_all <- ggplot() +
  geom_tile(data = entropy, aes(x = n, y = beta, fill = Dind_mean, colour = Dind_mean)) +
  theme_bw() +
  geom_hline(yintercept = 1, size = 0.3, linetype = "dashed") +
  geom_line(data = analytical_data, aes(x = n, y = beta), size = 0.6) +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20)), 
                     expand = c(0,0),
                     limits = c(2.5, 102.5)) +
  scale_y_continuous(breaks = seq(0.75, 1.25, 0.05), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste( "Interaction bias (", italic(beta), ")"))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 2)
gg_beta_all

ggsave(gg_beta_all, file = "output/ParameterSpace/Plots/BetaSweep_PositiveAll.png", height = 90, width = 45, units = "mm", dpi = 400)
ggsave(gg_beta_all, file = "output/ParameterSpace/Plots/svg/BetaSweep_PositiveAll.svg", height = 90, width = 45, units = "mm")

####################
# Plot: Beta sweep - negative social influence
####################
load("output/ParameterSpace/GroupSizeBetaSweep_Sigma0-Epsilon-0.1.Rdata")
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)

for (bias in c("Homophily", "Heterophily")) {
  # Filter
  if (bias == "Homophily") {
    entropy_filt <- entropy %>% 
      filter(beta >= 1)
  } else {
    entropy_filt <- entropy %>% 
      filter(beta <= 1)
  }
  # Graph
  gg_beta <- ggplot(data = entropy_filt, aes(x = n, y = beta, fill = Dind_mean, colour = Dind_mean)) +
    geom_tile() +
    theme_bw() +
    scale_x_continuous(breaks = c(5, seq(20, 100, 20)), 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0.75, 1.25, 0.05), 
                       expand = c(0,0)) +
    scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
    scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                           limits = c(0, 1)) +
    xlab(expression(paste("Group Size (", italic(n), ")"))) +
    ylab(expression(paste( "Interaction bias (", italic(beta), ")"))) +
    theme(axis.text = element_text(colour = "black", size = 6),
          axis.title = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.height = unit(5, "mm"),
          legend.key.width = unit(2, "mm"),
          legend.position = "none",
          axis.ticks = element_line(size = 0.3, color = "black"),
          panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
          aspect.ratio = 1)
  gg_beta
  # Save
  file_png <- paste0("output/ParameterSpace/Plots/BetaSweep_Negative", bias,  ".png")
  file_svg <- paste0("output/ParameterSpace/Plots/svg/BetaSweep_Negative", bias,  ".svg")
  ggsave(gg_beta, file = file_png, height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_beta, file = file_svg, height = 45, width = 45, units = "mm")
}

# Analytical results for group size above which full DOL should emerge (n*)
# Not accounting for double interact
analytical_data <- data.frame(beta = seq(1.001, 1.255, 0.0001), n = rep(NA, length(seq(1.001, 1.255, 0.0001))))
for (i in 1:nrow(analytical_data)) {
  analytical_data[i, 2] <- (2 * analytical_data[i, 1]) / (0.8 * (analytical_data[i, 1] - 1))
}
analytical_data <- analytical_data %>% 
  filter(n < 102.5)


# As one plot
gg_beta_all <- ggplot() +
  geom_tile(data = entropy, aes(x = n, y = beta, fill = Dind_mean, colour = Dind_mean)) +
  theme_bw() +
  geom_hline(yintercept = 1, size = 0.3, linetype = "dashed") +
  geom_line(data = analytical_data, aes(x = n, y = beta), size = 0.6) +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20)), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0.75, 1.25, 0.05), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste( "Interaction bias (", italic(beta), ")"))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 2)
gg_beta_all

ggsave(gg_beta_all, file = "output/ParameterSpace/Plots/BetaSweep_NegativeAll.png", height = 90, width = 45, units = "mm", dpi = 400)
ggsave(gg_beta_all, file = "output/ParameterSpace/Plots/svg/BetaSweep_NegativeAll.svg", height = 90, width = 45, units = "mm")


# Epsilon sweeps -----------------------------------------------------------------

####################
# Plot: Epsilon sweep - homophily
####################
load("output/ParameterSpace/GroupSizeEpsilonSweep_Sigma0-Beta1.1.Rdata")
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)

for (influence in c("Positive", "Negative")) {
  # Filter
  if (influence == "Positive") {
    entropy_filt <- entropy %>% 
      filter(epsilon >= 0)
  } else {
    entropy_filt <- entropy %>% 
      filter(epsilon <= 0)
  }
  # Graph
  gg_eps <- ggplot(data = entropy_filt, aes(x = n, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
    geom_tile() +
    theme_bw() +
    scale_x_continuous(breaks = c(5, seq(20, 100, 20)), 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(-1, 1, 0.1), 
                       expand = c(0,0)) +
    scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
    scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                           limits = c(0, 1)) +
    xlab(expression(paste("Group Size (", italic(n), ")"))) +
    ylab(expression(paste( "Social influence (", italic(epsilon), ")"))) +
    theme(axis.text = element_text(colour = "black", size = 6),
          axis.title = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.height = unit(5, "mm"),
          legend.key.width = unit(2, "mm"),
          legend.position = "none",
          axis.ticks = element_line(size = 0.3, color = "black"),
          panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
          aspect.ratio = 1)
  
  gg_eps
  # Save
  file_png <- paste0("output/ParameterSpace/Plots/EpsilonSweep_", influence,  "Homophily.png")
  file_svg <- paste0("output/ParameterSpace/Plots/svg/EpsilonSweep_", influence,  "Homophily.svg")
  ggsave(gg_eps, file = file_png, height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_eps, file = file_svg, height = 45, width = 45, units = "mm")
}

# As one plot
gg_eps_all <- ggplot(data = entropy, aes(x = n, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
  geom_tile() +
  geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20)), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste( "Social influence (", italic(epsilon), ")"))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 2)

gg_eps_all
ggsave(gg_eps_all, file = "output/ParameterSpace/Plots/EpsilonSweep_HomophilyAll.png", height = 90, width = 45, units = "mm", dpi = 400)
ggsave(gg_eps_all, file = "output/ParameterSpace/Plots/svg/EpsilonSweep_HomophilyAll.svg", height = 90, width = 45, units = "mm")


####################
# Plot: Epsilon sweep - heterophily
####################
load("output/ParameterSpace/GroupSizeEpsilonSweep_Sigma0-Beta0.9.Rdata")
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)

for (influence in c("Positive", "Negative")) {
  # Filter
  if (influence == "Positive") {
    entropy_filt <- entropy %>% 
      filter(epsilon >= 0)
  } else {
    entropy_filt <- entropy %>% 
      filter(epsilon <= 0)
  }
  # Graph
  gg_eps <- ggplot(data = entropy_filt, aes(x = n, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
    geom_tile() +
    theme_bw() +
    scale_x_continuous(breaks = c(5, seq(20, 100, 20)), 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(-1, 1, 0.1), 
                       expand = c(0,0)) +
    scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
    scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                           limits = c(0, 1)) +
    xlab(expression(paste("Group Size (", italic(n), ")"))) +
    ylab(expression(paste( "Social influence (", italic(epsilon), ")"))) +
    theme(axis.text = element_text(colour = "black", size = 6),
          axis.title = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.height = unit(5, "mm"),
          legend.key.width = unit(2, "mm"),
          legend.position = "none",
          axis.ticks = element_line(size = 0.3, color = "black"),
          panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
          aspect.ratio = 1)
 
  
  gg_eps
  # Save
  file_png <- paste0("output/ParameterSpace/Plots/EpsilonSweep_", influence,  "Heterophily.png")
  file_svg <- paste0("output/ParameterSpace/Plots/svg/EpsilonSweep_", influence,  "Heterophily.svg")
  ggsave(gg_eps, file = file_png, height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_eps, file = file_svg, height = 45, width = 45, units = "mm")
}

# As one plot
gg_eps_all <- ggplot(data = entropy, aes(x = n, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
  geom_tile() +
  geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(breaks = c(5, seq(20, 100, 20)), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.1), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste( "Social influence (", italic(epsilon), ")"))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 2)

gg_eps_all

ggsave(gg_eps_all, file = "output/ParameterSpace/Plots/EpsilonSweep_HeterophilyAll.png", height = 90, width = 45, units = "mm", dpi = 400)
ggsave(gg_eps_all, file = "output/ParameterSpace/Plots/svg/EpsilonSweep_HeterphilyAll.svg", height = 90, width = 45, units = "mm")



# Beta- Epsilon sweeps -----------------------------------------------------------------
 
####################
# Load
####################
load("output/ParameterSpace/EpsilonBetaSweep-n80.Rdata")
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)


####################
# Plot: Positive Homophily only
####################

entropy_filtered <- entropy %>% 
  filter(beta >= 1, 
         epsilon >= 0)

load("output/AnalyticalResults/EpsStar-BetaStar_Calc.Rdata")

# Graph
gg_betaeps <- ggplot() +
  geom_tile(data = entropy_filtered, aes(x = beta, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
  geom_line(data = eps_star, aes(x = beta, y = epsilon_all), 
            size = 0.3, 
            color = "grey40") +
  geom_vline(xintercept = unique(beta_star$beta), 
             size = 0.3, 
             color = "grey40") + 
  geom_hline(yintercept = 0, 
             size = 0.3, 
             color = "grey40") + 
  theme_bw() +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05), 
                     expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), 
                     expand = c(0,0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  xlab(expression(paste("Interaction Bias (", italic(beta), ")"))) +
  ylab(expression(paste( "Social influence (", italic(epsilon), ")"))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 1)


gg_betaeps
ggsave(gg_betaeps, file = "output/ParameterSpace/Plots/BetaEpsSweep_n80.png", height = 45, width = 45, units = "mm", dpi = 400)


ggsave(gg_betaeps, file = "output/ParameterSpace/Plots/svg/BetaEpsSweep_n80.svg", height = 45, width = 45, units = "mm")
ggsave(gg_betaeps, file = "output/ParameterSpace/Plots/svg/BetaEpsSweep_n80_large.svg", height = 75, width = 75, units = "mm")


####################
# Plot: All combos
####################

# Filter data


# Plot
gg_betaeps_all <- ggplot(data = entropy, aes(x = beta, y = epsilon, fill = Dind_mean, colour = Dind_mean)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0.75, 1.25, 0.05), 
                     expand = c(0,0),
                     labels = c("0.75", rep("", length(seq(0.75, 1, 0.05))-2), "1.00",rep("", length(seq(0.75, 1, 0.05))-2), "1.25")) +
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.1), 
                     expand = c(0,0),
                     labels = c("-0.6", rep("", length(seq(-0.6, 0, 0.1))-2), "0.0", rep("", length(seq(-0.6, 0, 0.1))-2), "0.6")) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  scale_colour_gradientn(colours = pal, name = "Behavioral\nspecialization",
                         limits = c(0, 1)) +
  xlab(expression(paste("Interaction Bias (", italic(beta), ")"))) +
  ylab(expression(paste( "Social influence (", italic(epsilon), ")"))) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
        aspect.ratio = 1) +
  geom_vline(xintercept = 1, size = 0.3, linetype = "dotted") +
  geom_vline(xintercept = 1.032258, size = 0.3) + #beta star
  geom_hline(yintercept = 0, size = 0.3, linetype = "dotted")
gg_betaeps_all

# ggsave(gg_betaeps_all, file = "output/ParameterSpace/Plots/BetaEpsSweep_n80_allcombos.png", height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_betaeps_all, file = "output/ParameterSpace/Plots/svg/BetaEpsSweep_n80_allcombos.svg", height = 75, width = 75, units = "mm")
