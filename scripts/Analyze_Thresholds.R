################################################################################
#
# Analyzing thresholds (at end point)
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)
library(viridis)
library(ggridges)

p <- 1 #prob of interact
runs <- c("Sigma0-Epsilon0.1-Beta1.1",
          "Sigma0.05-Epsilon0-Beta1.1")

####################
# Load and process data
####################
# Loop through group sizes
thresh_data <- lapply(runs, function(run) {
  # List files 
  files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh/", run, "/"), full.names = TRUE)
  # Loop through group sizes
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
  if (grepl(".*Epsilon0-.*", run)) {
    all_thresh$Model <- "Fixed thresholds"
  } else {
    all_thresh$Model <- "Socially-modulated thresholds"
  }
  return(all_thresh)
})

# bind and format
thresh_data <- as.data.frame(do.call('rbind', thresh_data))

####################
# Plot threshold distributions
####################
# Filter
filter_data <- thresh_data %>% 
  filter(n %in% seq(5, 50, 5))

# Plot
gg_threshvar <- ggplot(data = filter_data, 
                       aes(x = Thresh1, y = n, fill = Model, group = n)) +
  theme_invisible() +
  geom_density_ridges(size = 0.1, stat = "binline", bins = 100) +
  xlab(expression(paste("Threshold 1 value (", italic(theta[i1]), ")"))) +
  ylab(expression(paste("Group Size (", italic(n), ")"))) +
  scale_x_continuous(breaks = seq(0, 100, 25), 
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

gg_threshvar

ggsave(gg_threshvar, file = paste0("output/Thresholds/GroupSize_Thresholds.png"), width = 90, height = 60, units = "mm", dpi = 600)
ggsave(gg_threshvar, file = paste0("output/Thresholds/GroupSize_Thresholds.svg"), width = 90, height = 60, units = "mm", dpi = 600)

####################
# Plot frequency of threshold values within 1, 2 SD of mean
####################
# Calcualte values signifying mean, 1 SD, 2SD
thresh_dist_value <- all_thresh %>% 
  group_by(n) %>% 
  mutate(within_1SD = ( mean(Thresh1) - sd(Thresh1) < Thresh1 )  & Thresh1 < (mean(Thresh1) + sd(Thresh1))) %>% 
  mutate(within_2SD = ( mean(Thresh1) - 2*sd(Thresh1) < Thresh1 )  & Thresh1 < (mean(Thresh1) + 2*sd(Thresh1)) & within_1SD == FALSE) %>% 
  mutate(outside_2SD = within_1SD == FALSE & within_2SD == FALSE) %>% 
  summarise(within_1SD = sum(within_1SD),
            within_2SD = sum(within_2SD),
            outside_2SD = sum(outside_2SD)) %>% 
  mutate(within_1SD = within_1SD / (n*100),
         within_2SD = within_2SD / (n*100),
         outside_2SD = outside_2SD / (n*100))



####################
# Plot thresholds by replicate
####################
look <- all_thresh %>% 
  filter(n == 10) %>% 
  mutate(replicate = paste(sim, chunk, sep = "-"))

gg_thresh <- ggplot(look, aes(y = replicate, x = Thresh2,  color = replicate)) +
  geom_point() +
  scale_x_continuous(limits = c(39, 61)) +
  theme_classic() +
  theme(legend.position = "none")
gg_thresh
