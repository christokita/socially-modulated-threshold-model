################################################################################
#
# Comparing various specialization plots
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


####################
# Compare entropies
####################
# Load social
load("output/Sigma001-Eps001-Phi001-ConnectP01-Bias1.1.Rdata")

entropy <- unlist(groups_entropy, recursive = FALSE)
entropy <- do.call("rbind", entropy)  %>% 
  mutate(set = paste(n, replicate, sep = "-"))%>% 
  select(-Dsym, -Dtask) %>% 
  filter(n != 1) %>% 
  group_by(n) %>% 
  summarise(Mean = mean(Dind),
            SE = sd(Dind) / sqrt(length(Dind))) %>% 
  mutate(Model = "Social")

# Load non-social
load("output/Sigma001-Eps0-Phi0.Rdata")

entropy1 <- unlist(groups_entropy, recursive = FALSE)
entropy1 <- do.call("rbind", entropy1)  %>% 
  mutate(set = paste(n, replicate, sep = "-"))%>% 
  select(-Dsym, -Dtask) %>% 
  filter(n != 1) %>% 
  group_by(n) %>% 
  summarise(Mean = mean(Dind),
            SE = sd(Dind) / sqrt(length(Dind))) %>% 
  mutate(Model = "Fixed")

# Join
entropy <- rbind(entropy, entropy1)

####################
# Plot
####################
gg_entropy <- ggplot(data = entropy, aes(x = n, group = Model)) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE, color = Model), 
                width = 0.7) +
  geom_line(aes(y = Mean, color = Model, linetype = Model)) +
  geom_point(aes(y = Mean, color = Model), 
             size = 1.5,
             shape = 21) +
  theme_classic() +
  xlab("Group Size") +
  ylab("DOL Entropy") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = unique(entropy$n)) +
  scale_color_manual(values = c("black", "mediumseagreen")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(legend.position = "right", 
        legend.title = element_text(size = 7, face = "bold"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width= unit(0.4, "cm"),
        legend.margin =  margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, "cm"),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank(),
        strip.text = element_text(size = 8, face = "italic"),
        strip.background = element_rect(fill = NA, colour = NA),
        panel.spacing = unit(0.5, "cm"))

gg_entropy

ggsave(gg_entropy, file = "output/SpecializationPlots/Sigma001-Eps001-Phi001-ConnectP01-Bias1.1.png", height = 3, width = 3.5, units = "in", dpi = 800)

