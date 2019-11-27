################################################################################
#
# Analyzing behavior by looking at frequency of "specialists" and "generalists"
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)
library(viridis)
library(ggridges)


####################
# Load data
####################
load("output/Rdata/_ProcessedData/TaskDist/Sigma0-Epsilon0.1-Beta1.1.Rdata")
social_data <- compiled_data
social_data$Model <- "Socially-modulated thresholds"

load("output/Rdata/_ProcessedData/TaskDist/Sigma0.05-Epsilon0-Beta1.1.Rdata")
fixed_data <- compiled_data
fixed_data$Model <- "Fixed threhsolds"

behav_data <- rbind(social_data, fixed_data) %>% 
  mutate(Set = paste(sim, chunk, sep = "-")) %>% 
  mutate(Task_bias = Task1 - Task2,
         Task_ratio = Task1 - Task2 / (Task1 + Task2),
         Activity = Task1 + Task2)
rm(compiled_data, social_data, fixed_data)

####################
#  Look at raw behavioral distributions
####################
behav_analysis <- behav_data %>% 
  filter(n %in% c(10, 25, 50, 90))


# Plot behavioral distributions
gg_behav_dist <- ggplot(data = behav_analysis, aes(x = Task1, y = Task2)) +
  geom_point(size = 0.1, alpha = 0.01) +
  theme_ctokita() +
  facet_grid(Model~n)
gg_behav_dist

# Plot task bias vs activity level
gg_behav_act <- ggplot(data = behav_analysis, aes(x = Task_bias, y = Activity)) +
  geom_point(size = 0.1, alpha = 0.01) +
  theme_ctokita() +
  facet_grid(Model~n)
gg_behav_act

####################
# Behavior distribution by group size
####################
# Filter to 
filter_data <- behav_data %>% 
  filter(n %in% seq(5, 50, 5))

gg_behavvar<- ggplot(data = filter_data, 
                     aes(x = Task1, y = n, fill = Model, group = n)) +
  theme_invisible() +
  geom_density_ridges2(size = 0.1, stat = "binline", bins = 100) +
  xlab(expression(paste("Task 1 performance frequency (", italic(x[i1]), ")"))) +
  ylab(expression(paste("Group Size (", italic(n), ")"))) +
  scale_x_continuous(breaks = seq(-1, 1, 0.2), 
                     # limits = c(0, 1),
                     expand = c(0.03, 0)) +
  scale_y_continuous(breaks = c(5, seq(10, 50, 10)),
                     expand = c(0.03, 0)) +
  # scale_fill_viridis() +
  # scale_color_viridis() +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, color = "black"),
        axis.title.y = element_text(size = 7, color = "black", face = 'italic'),
        legend.position = "none",
        strip.text = element_text(size = 6, face = "bold"),
        axis.ticks = element_line(size = 0.2, color = "black")) +
  facet_grid(~Model)

gg_behavvar

ggsave(gg_behavvar, file = "output/Behavior/GroupSize_Behavior.png", width = 90, height = 60, units = "mm", dpi = 600)
ggsave(gg_behavvar, file = "output/Behavior/GroupSize_Behavior.svg", width = 90, height = 60, units = "mm", dpi = 600)


####################
#  Look at frequency of specialists and generalists
####################
# Categorize each individual as specialists or generalist
# 
# Def: 
# Generalist: >80% of activity on one task
#
behav_data <- behav_data %>% 
  mutate(Type = ifelse(abs(Task_ratio) > 0.8 & (Activity > 0), "Specialist", "Generalist"),
         Activity_type = ifelse(Activity > 0.5, "HighActivity", "LowActivity")) %>% 
  mutate(Type_refined = paste0(Type, "-", Activity_type))

behav_summary <- behav_data %>% 
  group_by(n, Model, Set, Type) %>% 
  summarise(Count = n()) %>% 
  mutate(Freq = Count / n) %>% 
  group_by(n, Model, Type) %>% 
  summarise(Mean_Count = mean(Count),
            SD_Count = sd(Count),
            Mean_Freq = mean(Freq),
            SD_Freq = sd(Freq))

# Plot
gg_type_freq <- ggplot(behav_summary, aes(x = n, y = Mean_Count, colour = Type)) +
  geom_line(size = 0.4) +
  geom_errorbar(aes(ymin = Mean_Count - SD_Count, ymax = Mean_Count + SD_Count), 
                width = 0,
                size = 0.3) +
  geom_point(size = 0.8) +
  scale_colour_manual(name = "Individual type", 
                      values = c("#a6cee3", "#1f78b4")) +
  ylab("Number of individuals") +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  scale_x_continuous(breaks = seq(20, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
  theme_ctokita() +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        legend.position = c(0.2, 0.75),
        legend.title = element_text(size = 6, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        strip.text = element_text(size = 6, face = "bold"),
        axis.ticks = element_line(size = 0.2, color = "black"),
        strip.background = element_rect(colour = NA, fill = "#DCDCDC")) +
  facet_grid(~Model) 

gg_type_freq
ggsave(gg_type_freq, file = "output/Behavior/GroupSize_BehaviorType.svg", width = 90, height = 45, units = "mm")

