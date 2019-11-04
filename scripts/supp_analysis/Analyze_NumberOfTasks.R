################################################################################
#
# Determine whether task order (i.e., order in which stimuli are encountered) affects simiulation results
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

############### m = 3, random stim encounter ############### 

####################
# Load data
####################
# Load data from runs with 3 tasks
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep_m3.Rdata")
data_m3 <- compiled_data %>% 
  mutate(TaskNumber = "Three") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

# Load data from runs with 2 tasks
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep.Rdata")
data_m2 <- compiled_data %>% 
  mutate(TaskNumber = "Two") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

rm(compiled_data)

# Bind together
tasknumber <- data_m3 %>% 
  rbind(data_m2)

####################
# Plot
####################
gg_tasknumber <- ggplot(data = tasknumber, aes(x = beta, color = TaskNumber)) +
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
  scale_color_manual(name = "Number of tasks", values = c("#ec7014", "#4d4d4d")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme_ctokita() +
  theme(axis.title.y = element_text(vjust = -1.5))
gg_tasknumber



############### m = 3, fixed and random stim encounter ############### 

####################
# Load data
####################
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


# Load data from runs with 2 tasks
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1_BetaSweep.Rdata")
data_m2 <- compiled_data %>% 
  mutate(TaskNumber = "Two") %>% 
  group_by(beta, TaskNumber) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

rm(compiled_data)

# Bind together
tasknumberorder <- data_m3 %>% 
  rbind(data_m2) %>% 
  rbind(data_m3order)

####################
# Plot
####################
gg_tasknumberorder <- ggplot(data = tasknumberorder, aes(x = beta, color = TaskNumber)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Interaction bias (", italic(beta), ")"))) +
  ylab(expression(paste("DOL (", italic(D[indiv]), ")"))) +
  scale_color_manual(name = "Number of tasks", values = c("#ec7014", "#fec44f", "#4d4d4d")) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5)) +
  theme_ctokita() +
  theme(axis.title.y = element_text(vjust = -1.5))
gg_tasknumberorder
