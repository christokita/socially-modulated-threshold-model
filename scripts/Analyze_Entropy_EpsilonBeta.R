################################################################################
#
# Comparing various specialization plots
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


############### Epsilon-Beta plots  ###############
rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

####################
# Load data
####################
load("output/ParameterSpace/EpsilonBetaSweep-n60.Rdata")

# Filter to parameter combos of interest
epsilons_of_interest <- c(0, 0.1, 0.3, 0.5)
betas_of_interest <- c(1, 1.05, 1.1, 1.2)

betas <- entropy %>% 
  filter(epsilon %in% epsilons_of_interest) %>% 
  mutate(Mean = Dind_mean,
         SD = Dind_SD,
         epsilon = as.factor(epsilon))

epsilons <- entropy %>% 
  filter(beta %in% betas_of_interest) %>% 
  mutate(Mean = Dind_mean,
         SD = Dind_SD,
         beta = as.factor(beta))

####################
# Plot betas
####################
# Betas
pal <- brewer.pal(5, "Greens")[2:5]

gg_entropy_betas <- ggplot(data = betas, aes(x = beta, colour = epsilon)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_color_manual(values = pal, 
                     labels = c("0.0", "0.1", "0.3", "0.5"),
                     name = expression("Social\ninfluence "(epsilon))) +
  theme_ctokita() +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = "none",
        legend.title = element_text(size = 6, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_entropy_betas

ggsave(gg_entropy_betas, file = "output/SpecializationPlots/Beta-Eps-n60-DiffEpsValues.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_entropy_betas, file = "output/SpecializationPlots/Beta-Eps-n60-DiffEpsValues.svg", 
       heigh = 45, width = 45, units = "mm")

# resize
gg_entropy_betas <- gg_entropy_betas + 
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme(aspect.ratio = NULL,
        axis.title.y = element_blank())

ggsave(gg_entropy_betas, file = "output/SpecializationPlots/Beta-Eps-n60-DiffEpsValues.svg", 
       width = 45, height = 25, units = "mm")

####################
# Plot epsilons
####################
# Epsilon
pal <- brewer.pal(5, "Greens")[2:5]

gg_entropy_eps <- ggplot(data = epsilons, aes(x = epsilon, colour = beta)) +
  geom_errorbar(aes(ymin = ifelse((Mean - SD) > 0 , Mean - SD, 0), ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Social influence (", italic(epsilon), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  scale_x_continuous(breaks = seq(0, 0.6, 0.1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = pal, 
                     labels = as.character(betas_of_interest),
                     name = expression("Interaction\nbias "(beta))) +
  theme_ctokita() +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = "right",
        legend.title = element_text(size = 6, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.ticks = element_line(size = 0.3),
        axis.line = element_line(size = 0.3),
        aspect.ratio = 1)
gg_entropy_eps

ggsave(gg_entropy_eps, file = "output/SpecializationPlots/Beta-Eps-n60-DiffBetaValues.png", 
       height = 45, width = 45, units = "mm", dpi = 400)
ggsave(gg_entropy_eps, file = "output/SpecializationPlots/Beta-Eps-n60-DiffBetaValues.svg", 
       width = 45, height = 45, units = "mm")

# resize
gg_entropy_eps <- gg_entropy_eps + 
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme(aspect.ratio = NULL,
        axis.title.y = element_blank())

ggsave(gg_entropy_eps, file = "output/SpecializationPlots/Beta-Eps-n60-DiffBetaValues.svg", 
       width = 45, height = 25, units = "mm")
