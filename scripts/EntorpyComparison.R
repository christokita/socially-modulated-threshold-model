################################################################################
#
# Comparing various specialization plots
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


filename <- "Sigma001-Eps001-Phi001-ConnectP01-Bias1.1"

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

ggsave(gg_entropy, file = paste0("output/SpecializationPlots/", filename, ".png"), height = 3, width = 3.5, units = "in", dpi = 800)


####################
# Create inset for task rank correlation
####################
# Load social
load("output/Sigma001-Eps001-Phi001-ConnectP01-Bias1.1.Rdata")

taskCorrTot <- do.call("rbind", groups_taskCorr)
taskCorrTot <- taskCorrTot %>% 
  mutate(TaskMean = (Task1 + Task2) / 2) %>% 
  mutate(set = paste(n, replicate, sep = "-"))%>% 
  select(-Task1, -Task2) %>% 
  group_by(n) %>% 
  summarise(Mean = mean(TaskMean),
            SE = sd(TaskMean) / sqrt(length(TaskMean))) %>% 
  mutate(Model = "Social")

# Load non-social
load("output/Sigma001-Eps0-Phi0.Rdata")

taskCorrTot1 <- do.call("rbind", groups_taskCorr)
taskCorrTot1 <- taskCorrTot1 %>% 
  mutate(TaskMean = (Task1 + Task2) / 2) %>% 
  mutate(set = paste(n, replicate, sep = "-"))%>% 
  select(-Task1, -Task2) %>% 
  group_by(n) %>% 
  summarise(Mean = mean(TaskMean),
            SE = sd(TaskMean) / sqrt(length(TaskMean))) %>% 
  mutate(Model = "Fixed")

# Join
taskCorrTot <- rbind(taskCorrTot, taskCorrTot1)


gg_corr <- ggplot(data = taskCorrTot, aes(x = n, group = Model)) +
  geom_line(aes(y = Mean, color = Model, linetype = Model)) +
  theme_classic() +
  xlab("Group Size") +
  ylab("Rank Correlation") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = unique(entropy$n)) +
  scale_color_manual(values = c("black", "mediumseagreen")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme(legend.position = "none", 
        legend.title = element_text(size = 7, face = "bold"),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width= unit(0.4, "cm"),
        legend.margin =  margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, "cm"),
        legend.text = element_text(size = 10),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 8, face = "italic"),
        strip.background = element_rect(fill = NA, colour = NA),
        panel.spacing = unit(0.5, "cm"))

gg_corr

ggsave(gg_corr, file = paste0("output/SpecializationPlots/", filename, "_TaskCorr.png"), height = 1, width = 1, units = "in", dpi = 800)

