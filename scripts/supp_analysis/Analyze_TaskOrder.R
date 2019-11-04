################################################################################
#
# Determine whether task order (i.e., order in which stimuli are encountered) affects simiulation results
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


####################
# Load data
####################
# Load data from runs with fixed order
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep_stimorder.Rdata")
data_fixedorder <- compiled_data %>% 
  mutate(Order = "Fixed") %>% 
  group_by(beta, Order) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

# Load data from runs with random order (i.e., runs normally used in analysis in main text)
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep.Rdata")
data_randorder <- compiled_data %>% 
  mutate(Order = "Random") %>% 
  group_by(beta, Order) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

rm(compiled_data)

# Bind together
stimorder <- data_randorder %>% 
  rbind(data_fixedorder)

####################
# Plot
####################
gg_stimorder <- ggplot(data = stimorder, aes(x = beta, color = Order)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3,
                position = position_dodge(width = 0.01)) +
  geom_point(aes(y = Mean),
             size = 0.8,
             position = position_dodge(width = 0.01)) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_color_manual(name = "Order of stimulus\nencounter", values = c("#ec7014", "#4d4d4d")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme_ctokita() +
  theme(axis.title.y = element_text(vjust = -1.5))
gg_stimorder
