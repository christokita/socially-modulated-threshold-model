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
# Set plotting theme for this script
####################
plot_settings <- function() {
  theme(axis.title.y = element_text(vjust = -1.5),
                        legend.position = c(0.85, 0.2),
                        legend.key.height = unit(2, "mm"),
                        legend.background = element_blank())
}

m2color <- "#4d4d4d"
m3color <- "#fec44f"
m5color <- "#ec7014"
subplt_height <- 21
subplt_legendpos <- c(0.75, 0.4)


####################
# Load data
####################
# Load data from runs with 2 tasks
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep.Rdata")
data_m2 <- compiled_data %>% 
  mutate(TaskNumber = "Two") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

# Load data from runs with 2 tasks, fixed order
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep_stimorder.Rdata")
data_m2order <- compiled_data %>% 
  mutate(TaskNumber = "Two, Fixed Order") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

# Load data from runs with 3 tasks
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep_m3.Rdata")
data_m3 <- compiled_data %>% 
  mutate(TaskNumber = "Three") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

# Load data from runs with 3 tasks, fixed order
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep_m3stimorder.Rdata")
data_m3order <- compiled_data %>% 
  mutate(TaskNumber = "Three, Fixed Order") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

# Load data from runs with 5 tasks
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep_m5.Rdata")
data_m5 <- compiled_data %>% 
  mutate(TaskNumber = "Five") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

# Load data from runs with 5 tasks, fixed order
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep_m5stimorder.Rdata")
data_m5order <- compiled_data %>% 
  mutate(TaskNumber = "Five, Fixed Order") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

rm(compiled_data)

# Bind together
tasknumberorder <- data_m2 %>% 
  rbind(data_m2order) %>% 
  rbind(data_m3) %>% 
  rbind(data_m3order) %>% 
  rbind(data_m5) %>% 
  rbind(data_m5order) %>% 
  mutate(TaskNumber = factor(TaskNumber, levels = c("Two", "Two, Fixed Order",
                                                    "Three", "Three, Fixed Order", 
                                                    "Five", "Five, Fixed Order")))


####################
# Plot, task number (all random stim encounter)
####################
# Filter data
tasknumber <- tasknumberorder %>% 
  filter(TaskNumber %in% c("Two", "Three", "Five"))

# Plot and save
gg_tasknumber <- ggplot(data = tasknumber, aes(x = beta, color = TaskNumber)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_color_manual(name = "", 
                     values = c(m2color, m3color, m5color), 
                     labels = c("m = 2", "m = 3", "m = 5")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme_ctokita() +
  plot_settings() +
  theme(aspect.ratio = 1)
ggsave(gg_tasknumber, filename = "output/TaskNumber/DOL_tasknumber.png", width = 45, height = 45, units = "mm", dpi = 400)
ggsave(gg_tasknumber, filename = "output/TaskNumber/DOL_tasknumber.svg", width = 45, height = 45, units = "mm")

####################
# m = 2, task encounter order
####################
# Filter data
m2 <- tasknumberorder %>% 
  filter(TaskNumber %in% c("Two", "Two, Fixed Order"))

# Plot and save
gg_m2taskorder <- ggplot(data = m2, aes(x = beta, color = TaskNumber, shape = TaskNumber)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_color_manual(name = "", 
                     values = c(m2color, m2color), 
                     labels = c("m = 2", "m = 2, fixed order")) +
  scale_shape_manual(name = "", 
                     values = c(16, 1), 
                     labels = c("m = 2", "m = 2, fixed order")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 1)) +
  theme_ctokita() +
  plot_settings() +
  theme(axis.title.y = element_blank(),
        legend.position = subplt_legendpos)
ggsave(gg_m2taskorder, filename = "output/TaskNumber/DOL_m2taskorder.svg", width = 45, height = subplt_height, units = "mm")

# larger plot for revision letter
gg_m2revisletter <- ggplot(data = m2, aes(x = beta, color = TaskNumber, shape = TaskNumber)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3,
                position = position_dodge(width = 0.015)) +
  geom_point(aes(y = Mean),
             size = 0.8,
             position = position_dodge(width = 0.015)) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_color_manual(name = "", 
                     values = c(m2color, m2color), 
                     labels = c("m = 2", "m = 2, fixed order")) +
  scale_shape_manual(name = "", 
                     values = c(16, 1), 
                     labels = c("m = 2", "m = 2, fixed order")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 1)) +
  theme_ctokita() +
  plot_settings() +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.72, 0.2))
ggsave(gg_m2revisletter, filename = "output/TaskNumber/DOL_m2taskorder.png", width = 45, height = 45, units = "mm", dpi = 400)


####################
# m = 3, task encounter order
####################
# Filter data
m3 <- tasknumberorder %>% 
  filter(TaskNumber %in% c("Three", "Three, Fixed Order"))

# Plot and save
gg_m3taskorder <- ggplot(data = m3, aes(x = beta, color = TaskNumber, shape = TaskNumber)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_color_manual(name = "", 
                     values = c(m3color, m3color), 
                     labels = c("m = 3", "m = 3, fixed order")) +
  scale_shape_manual(name = "", 
                     values = c(16, 1), 
                     labels = c("m = 3", "m = 3, fixed order")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 1)) +
  theme_ctokita() +
  plot_settings() +
  theme(axis.title.y = element_blank(),
        legend.position = subplt_legendpos)
ggsave(gg_m3taskorder, filename = "output/TaskNumber/DOL_m3taskorder.svg", width = 45, height = subplt_height, units = "mm")

####################
# m = 5, task encounter order
####################
# Filter data
m5 <- tasknumberorder %>% 
  filter(TaskNumber %in% c("Five", "Five, Fixed Order"))

# Plot and save
gg_m5taskorder <- ggplot(data = m5, aes(x = beta, color = TaskNumber, shape = TaskNumber)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_color_manual(name = "", 
                     values = c(m5color, m5color), 
                     labels = c("m = 5", "m = 5, fixed order")) +
  scale_shape_manual(name = "", 
                     values = c(16, 1), 
                     labels = c("m = 5", "m = 5, fixed order")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 1)) +
  theme_ctokita() +
  plot_settings() +
  theme(axis.title.y = element_blank(),
        legend.position = subplt_legendpos)
ggsave(gg_m5taskorder, filename = "output/TaskNumber/DOL_m5taskorder.svg", width = 45, height = subplt_height, units = "mm")

