################################################################################
#
# Testing whether the threshold limit is causing artifacts for high social influence
#
################################################################################
#
# Here we are testing wehther the decrease in DOL at high levels of social influence
# is due to the upper limit on threhsold values. In this simulation, we put the upper bound
# for threshold at 1000, instead of the normal 100.

rm(list = ls())
source("scripts/util/__Util__MASTER.R")

# ------------------------------ DOL ------------------------------
####################
# Load and process data
####################
load('output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.4-Beta1.1-HighThreshLimit.Rdata')
high_thresh <- compiled_data
high_thresh$Model <- "high_thresh"

load('output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.4-Beta1.1.Rdata')
normal_thresh <- compiled_data
normal_thresh$Model <- "normal_thresh"

entropy_data <- rbind(high_thresh, normal_thresh)
entropy_data <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            SD = sd(Dind))

####################
# Plot entropy plots
####################
gg_comp <- ggplot(entropy_data, aes(x = n, y = Mean, group = Model, color = Model)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0,
                size = 0.3) +
  geom_line(aes(y = Mean),
            size = 0.4) +
  geom_point(aes(y = Mean),
             size = 0.8) +
  theme_classic() +
  xlab(expression(paste("Group Size (", italic(n), ")"))) +
  ylab(expression(paste("Division of labor (", italic(D[indiv]), ")"))) +
  ggtitle(expression(paste(italic(epsilon), "= 0.4, ", italic(beta), "= 1.1"))) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme(title = element_text(size = 6),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, face = "italic"),
        legend.position = c(0.8, 0.2),
        legend.title = element_text(size = 6, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)
gg_comp
ggsave(gg_comp, file = "output/SpecializationPlots/ThresholdLimitComparison.png", height = 50, width = 50, units = "mm")
  

# ------------------------------ Thresholds ------------------------------
####################
# Load and process data
####################
files <- list.files('output/Rdata/_ProcessedData/Thresh/Sigma0-Epsilon0.4-Beta1.1-HighThreshLimit/', full.names = T)
group_thresh <- lapply(files, function(file) {
  
  # Load
  load(file)
  
  # Bind
  thresh_data <- as.data.frame(do.call("rbind", listed_data))
  
  # # Summarise
  # thresh_sum <- thresh_data %>% 
  #   mutate(run = paste0(sim, "-", chunk)) %>% 
  #   group_by(n, run) %>% 
  #   summarise(Thresh1SD = sd(Thresh1),
  #             Thresh2SD = sd(Thresh2),
  #             ThreshSD = (sd(Thresh1)  + sd(Thresh2)) / 2) %>% 
  #   select(n, ThreshSD)
  return(thresh_data)
})
# Bind 
all_thresh <- do.call('rbind', group_thresh)


####################
# Plot distirubtion
####################
library(ggridges)
library(viridis)
# Plot
gg_threshvar <- ggplot(data = all_thresh, 
                       aes(x = Thresh1, y = n, fill = n, group = n)) +
  theme_invisible() +
  geom_density_ridges(size = 0.1, stat = "binline", bins = 100) +
  xlab(expression(paste("Task 1 threshold (", italic(theta[i1]), ")"))) +
  ylab(expression(paste("Group Size (", italic(n), ")"))) +
  # scale_x_continuous(breaks = seq(0, 100, 25), 
  #                    # limits = c(0, 1),
  #                    expand = c(0.03, 0)) +
  scale_y_continuous(breaks = c(5, seq(10, 50, 10)),
                     expand = c(0.03, 0)) +
  scale_fill_viridis() +
  scale_color_viridis() +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, color = "black"),
        axis.title.y = element_text(size = 7, color = "black", face = 'italic'),
        legend.position = "none",
        strip.text = element_text(size = 6, face = "bold"),
        axis.ticks = element_line(size = 0.2, color = "black"))

gg_threshvar
