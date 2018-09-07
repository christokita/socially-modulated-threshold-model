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
# Plot: Beta sweep - positive social influence
####################
load("output/ParameterSpace/GroupSizeBetaSweep_Sigma0-Epsilon0.1.Rdata")
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
    scale_x_continuous(breaks = seq(0, 100, 20), 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0.75, 1.25, 0.05), 
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
  # Save
  file_png <- paste0("output/ParameterSpace/Plots/BetaSweep_Positive", bias,  ".png")
  file_svg <- paste0("output/ParameterSpace/Plots/svg/BetaSweep_Positive", bias,  ".svg")
  ggsave(gg_beta, file = file_png, height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_beta, file = file_svg, height = 45, width = 45, units = "mm")
}

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
  # Save
  file_png <- paste0("output/ParameterSpace/Plots/BetaSweep_Negative", bias,  ".png")
  file_svg <- paste0("output/ParameterSpace/Plots/svg/BetaSweep_Negative", bias,  ".svg")
  ggsave(gg_beta, file = file_png, height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_beta, file = file_svg, height = 45, width = 45, units = "mm")
}


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
  # Save
  file_png <- paste0("output/ParameterSpace/Plots/EpsilonSweep_", influence,  "Homophily.png")
  file_svg <- paste0("output/ParameterSpace/Plots/svg/EpsilonSweep_", influence,  "Homophily.svg")
  ggsave(gg_eps, file = file_png, height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_eps, file = file_svg, height = 45, width = 45, units = "mm")
}

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
  # Save
  file_png <- paste0("output/ParameterSpace/Plots/EpsilonSweep_", influence,  "Heterophily.png")
  file_svg <- paste0("output/ParameterSpace/Plots/svg/EpsilonSweep_", influence,  "Heterophily.svg")
  ggsave(gg_eps, file = file_png, height = 45, width = 45, units = "mm", dpi = 400)
  ggsave(gg_eps, file = file_svg, height = 45, width = 45, units = "mm")
}



