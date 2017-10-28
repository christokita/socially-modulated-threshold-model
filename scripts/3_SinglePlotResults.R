################################################################################
#
# Plot ensemble model outputs in single plot 
#
################################################################################
rm(list = ls())
source("scripts/__Util__MASTER.R")
source("scripts/3A_PrepPlotExperimentData.R")

# Set variable  
filename <- "Fixed_Delta06Sigma01Eta7"

# Palette without single individuals
#palette <- c("#F00924", "#F7A329", "#FDD545", "#027C2C", "#1D10F9", "#4C0E78", "#bdbdbd", "#525252")

# Palette without single individuals
palette <- c("#83343E", "#F00924", "#F7A329", "#FDD545", "#027C2C", "#1D10F9", "#4C0E78", "#bdbdbd", "#525252")

# Model vs Data Palette
compPalette <- c("indianred2", "black")

####################
# Final task distributions
####################
# Bind together
taskDist <- unlist(groups_taskDist, recursive = FALSE)
taskDistTot <- do.call("rbind", taskDist)

# Manipulate
taskDistTot <- taskDistTot %>% 
  mutate(set = paste0(n, "-", replicate)) %>% 
  mutate(set = factor(set, 
                      levels = mixedsort(unique(set))),
         n = as.factor(n))

taskSum <- taskDistTot %>% 
  group_by(n) %>% 
  summarise(taskMean1 = mean(Task1),
            taskMean2 = mean(Task2))

# Plot
plot_TaskMat <- as.data.frame(taskDistTot)
gg_dist <- ggplot(data = plot_TaskMat, aes(y = Task1, x = set)) +
  geom_point(aes(colour = n), size = 0.3) +
  theme_classic() +
  labs(x = "Group size",
       y = "Task 1 frequency") +
  scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme(axis.text.x = element_text(size = 0),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 10, margin = margin(0, 0, 0, 0)),
        axis.title.x = element_text(size = 10, margin = margin(12, 0, 0, 0)),
        legend.position = "none") 



####################
# Task variance by group size
####################
# Prep
taskVarMean <- taskDistTot %>% 
  mutate(n = as.character(n)) %>% 
  mutate(n = as.numeric(n)) %>% 
  group_by(n, replicate) %>% 
  summarise(SD1 = sd(Task1),
            SD2 = sd(Task2),
            Mean = mean(Task1)) %>% 
  mutate(Source = "Model",
         SD = (SD1 + SD2) / 2)
taskVarMean$SD[is.na(taskVarMean$SD)] <- 0 #fix for single individuals

taskVarMean <- rbind(taskVarMean, yukoDataSummary)


# Calculate means and SE
taskVarMeans <- taskVarMean %>% 
  group_by(n, Source) %>% 
  summarise(MeanMean = mean(Mean),
            MeanSE = sd(Mean) / sqrt(length(Mean)),
            SDMean = mean(SD),
            SDSE = sd(SD) / sqrt(length(SD)))

# Get mean at group size one and normalize
expSizeOne <- taskVarMeans$MeanMean[taskVarMeans$n == 1 & taskVarMeans$Source == "Experiment"]
modSizeOne <- taskVarMeans$MeanMean[taskVarMeans$n == 1 & taskVarMeans$Source == "Model"]

expSizeSixteen <- taskVarMeans$SDMean[taskVarMeans$n == 16 & taskVarMeans$Source == "Experiment"]
modSizeSixteen <- taskVarMeans$SDMean[taskVarMeans$n == 16 & taskVarMeans$Source == "Model"]

# Normalize Mean Values
taskVarMeans$NormMean <- NA
taskVarMeans$NormMeanSE <- NA
taskVarMeans$NormMean[taskVarMeans$Source == "Experiment"] <- taskVarMeans$MeanMean[taskVarMeans$Source == "Experiment"] / expSizeOne
taskVarMeans$NormMean[taskVarMeans$Source == "Model"] <- taskVarMeans$MeanMean[taskVarMeans$Source == "Model"] / modSizeOne
taskVarMeans$NormMeanSE[taskVarMeans$Source == "Experiment"] <- taskVarMeans$MeanSE[taskVarMeans$Source == "Experiment"] / expSizeOne
taskVarMeans$NormMeanSE[taskVarMeans$Source == "Model"] <- taskVarMeans$MeanSE[taskVarMeans$Source == "Model"] / modSizeOne

taskVarMean$NormMean <- NA
taskVarMean$NormMean[taskVarMean$Source == "Experiment"] <- taskVarMean$Mean[taskVarMean$Source == "Experiment"] / expSizeOne
taskVarMean$NormMean[taskVarMean$Source == "Model"] <- taskVarMean$Mean[taskVarMean$Source == "Model"] / modSizeOne

# Normalize behavioral variation values
taskVarMeans$NormVarMean <- NA
taskVarMeans$NormVarMeanSE <- NA
taskVarMeans$NormVarMean[taskVarMeans$Source == "Experiment"] <- taskVarMeans$SDMean[taskVarMeans$Source == "Experiment"] / expSizeSixteen
taskVarMeans$NormVarMean[taskVarMeans$Source == "Model"] <- taskVarMeans$SDMean[taskVarMeans$Source == "Model"] / modSizeSixteen
taskVarMeans$NormVarMeanSE[taskVarMeans$Source == "Experiment"] <- taskVarMeans$SDSE[taskVarMeans$Source == "Experiment"] / expSizeSixteen
taskVarMeans$NormVarMeanSE[taskVarMeans$Source == "Model"] <- taskVarMeans$SDSE[taskVarMeans$Source == "Model"] / modSizeSixteen

taskVarMean$NormVarMean <- NA
taskVarMean$NormVarMean[taskVarMean$Source == "Experiment"] <- taskVarMean$SD[taskVarMean$Source == "Experiment"] / expSizeSixteen
taskVarMean$NormVarMean[taskVarMean$Source == "Model"] <- taskVarMean$SD[taskVarMean$Source == "Model"] / modSizeSixteen


# Plot variance and mean by group size
# gg_var <- ggplot() +
#   geom_point(data = taskVarMean, 
#              aes(x = n, y = SD, colour = Source),
#              size = 0.5,
#              alpha = 0.4,
#              position = position_dodge(width = 1)) +
#   theme_classic() +
#   xlab("Group size") +
#   ylab("Behavioral Variation (SD)") +
#   scale_x_continuous(breaks = unique(taskVarMean$n)) +
#   scale_y_continuous(breaks = seq(0, 1, 0.025)) +
#   # Mean and SE portion of plot
#   geom_errorbar(data = taskVarMeans, 
#                 aes(x = n, ymin = SDMean - SDSE, ymax = SDMean + SDSE, colour = Source, width = 1.5),
#                 position = position_dodge(width = 1)) +
#   geom_point(data = taskVarMeans, 
#              aes(x = n, y = SDMean, colour = Source),
#              size = 2,
#              position = position_dodge(width = 1)) +
#   geom_line(data = taskVarMeans,
#             aes(x = n, y = SDMean, colour = Source),
#             position = position_dodge(width = 1)) +
#   scale_fill_manual(values = compPalette) +
#   scale_colour_manual(values = compPalette) +
#   theme(legend.position = "none",
#         axis.text = element_text(size = 8),
#         axis.title = element_text(size = 10))

gg_varNorm <- ggplot() +
  geom_hline(data = taskVarMean, 
             aes(yintercept = 1),
             colour = "grey30") +
  geom_point(data = taskVarMean, 
             aes(x = n, y = NormVarMean, colour = Source),
             size = 0.9,
             alpha = 0.4,
             position = position_dodge(width = 1),
             stroke = 0) +
  theme_classic() +
  xlab("Group size") +
  ylab("Relative behavioral variation") +
  scale_x_continuous(breaks = unique(taskVarMean$n)) +
  scale_y_continuous(breaks = seq(0, 3, 0.5)) +
  # Mean and SE portion of plot
  geom_errorbar(data = taskVarMeans, 
                aes(x = n, ymin = NormVarMean - NormVarMeanSE, ymax = NormVarMean + NormVarMeanSE, colour = Source),
                width = 1.5,
                position = position_dodge(width = 1)) +
  geom_point(data = taskVarMeans, 
             aes(x = n, y = NormVarMean, colour = Source),
             size = 2,
             position = position_dodge(width = 1)) +
  geom_line(data = taskVarMeans,
            aes(x = n, y = NormVarMean, colour = Source),
            position = position_dodge(width = 1)) +
  scale_fill_manual(values = compPalette) +
  scale_colour_manual(values = compPalette) +
  theme(legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10, margin = margin(0, 0, 0, 0)))


gg_mean <- ggplot() +
  geom_hline(data = taskVarMean, 
             aes(yintercept = 1),
             colour = "grey30") +
  geom_point(data = taskVarMean,
             aes(x = n, y = NormMean, colour = Source),
             size = 0.9,
             alpha = 0.4,
             position = position_dodge(width = 1),
             stroke = 0) +
  theme_classic() +
  xlab("Group size") +
  ylab("Relative task 1 frequency") +
  scale_x_continuous(breaks = unique(taskVarMean$n)) +
  scale_y_continuous(breaks = seq(0, 1.5, 0.05)) +
  # Mean and SE portion of plot
  geom_errorbar(data = taskVarMeans, 
                aes(x = n, ymin = NormMean - NormMeanSE, ymax = NormMean + NormMeanSE, colour = Source),
                width = 1.5,
                position = position_dodge(width = 1)) +
  geom_point(data = taskVarMeans, 
             aes(x = n, y = NormMean, colour = Source),
             size = 2,
             position = position_dodge(width = 1)) +
  geom_line(data = taskVarMeans,
            aes(x = n, y = NormMean, colour = Source),
            position = position_dodge(width = 1)) +
  scale_fill_manual(values = compPalette) +
  scale_colour_manual(values = compPalette) +
  theme(legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10, margin = margin(0, 0, 0, 0)))

####################
# Task Rank Correlation
####################
# Unlist
taskCorrTot <- do.call("rbind", groups_taskCorr)
taskCorrTot <- taskCorrTot %>% 
  mutate(TaskMean = (Task1 + Task2) / 2)

# Manipulate and bind with Yuko data
taskCorrTot <- taskCorrTot %>% 
  mutate(Source = "Model") %>% 
  select(n, TaskMean, Source) %>% 
  rbind(yukoCorr) %>% 
  mutate(Source = as.factor(Source))

# Calculate means and SE
taskCorrMeans <- taskCorrTot %>% 
  group_by(Source, n) %>% 
  summarise(SpecMean = mean(TaskMean),
            SpecSE = sd(TaskMean) / sqrt(length(TaskMean)),
            SpecCI = 1.96 * SpecSE)

# Plot
gg_corr <- ggplot() +
  geom_hline(data = taskCorrTot, 
             aes(yintercept = 0),
             colour = "grey30") +
  geom_point(data = taskCorrTot, 
             aes(x = n, y = TaskMean, fill = Source, colour = Source), 
             size = 0.9, 
             position = position_dodge(width = 1),
             alpha = 0.4,
             stroke = 0) +
  theme_classic() +
  labs(x = "Group size",
       y = "Specialization") +
  scale_x_continuous(breaks = unique(taskCorrTot$n)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
  scale_fill_manual(values = compPalette) +
  scale_colour_manual(values = compPalette) +
  theme(legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10, margin = margin(0, 0, 0, 0))) +
  # Mean and SE portion of plot
  geom_errorbar(data = taskCorrMeans, 
                aes(x = n, ymin = SpecMean - SpecSE, ymax = SpecMean + SpecSE, colour = Source),
                width = 1.5,
                position = position_dodge(width = 1)) +
  geom_point(data = taskCorrMeans, 
             aes(x = n, y = SpecMean, colour = Source, fill = Source),
             position = position_dodge(width = 1),
             size = 2) +
  geom_line(data = taskCorrMeans,
            aes(x = n, y = SpecMean,  colour = Source),
            position = position_dodge(width = 1))


####################
# Plot all
####################


# MultiPlot
png(filename = paste0("output/_ComprehnsivePlots/", filename, ".png"), width = 4, height = 4, units = "in", res = 800)
multiplot(gg_dist, gg_mean, gg_corr, gg_varNorm, cols = 2)  
dev.off()



####################
# Save all
####################
# save(groups_entropy, groups_stim, groups_taskCorr, groups_taskDist, groups_taskStep, groups_taskTally, taskCorrTot, 
#      file = paste0("output/__RData/", filename, ".Rdata"))

