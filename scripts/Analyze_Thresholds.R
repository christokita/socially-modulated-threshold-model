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
run <- "Sigma0-Epsilon0.1-Beta1.1"

####################
# Load and process data
####################
# Load social networks
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

# Format for plotting


# Plot
pal <- brewer_pal("seq", "RdPu")
pal <- pal(9)

gg_threshvar <- ggplot(data = all_thresh, 
                       aes(x = Thresh1, y = n, fill = n, group = n)) +
  theme_invisible() +
  geom_density_ridges(size = 0.2) +
  xlab("Threshold value") +
  scale_x_continuous(breaks = seq(0, 100, 25), 
                     expand = c(0.03, 0)) +
  scale_y_continuous(breaks = c(5, seq(25, 100, 25)), 
                     expand = c(0.03, 0)) +
  scale_fill_viridis() +
  scale_color_viridis() +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7, color = "black"),
        legend.position = "none",
        axis.ticks = element_line(size = 0.2, color = "black"))

gg_threshvar

ggsave(gg_threshvar, file = "output/Thresholds/GroupSizeThreshold.png", width = 90, height = 60, units = "mm", dpi = 600)
ggsave(gg_threshvar, file = "output/Thresholds/GroupSizeThreshold.svg", width = 90, height = 60, units = "mm", dpi = 600)


ggsave(gg_threshvar, file = "output/Thresholds/GroupSizeThreshold_square.png", width = 70, height = 70, units = "mm", dpi = 600)

