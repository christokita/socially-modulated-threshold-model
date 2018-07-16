################################################################################
#
# Analyzing thresholds over time
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)
library(viridis)
library(rje)

p <- 1 #prob of interact
run <- "Sigma0-Epsilon0.1-Beta1.1"

####################
# Load and process data
####################
# Load social networks
files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh1Time/", run, "/"), full.names = TRUE)

# Get group sizes
group_sizes <- unique(gsub(".*/([0-9]+)-[0-9]+\\.Rdata", "\\1", files, perl = TRUE))

select_sizes <- group_sizes[group_sizes %in% c("020", "025", 
                                              "030", "035",
                                              "050")]

# Loop through group sizes
group_plots <- lapply(select_sizes, function(size) {
  
  # Get subset of files that are for this group size
  group_files <- files[grep(paste0(".*/", size, "-[0-9]+\\.Rdata"), files)]
  
  # Loop through group files
  group_data <- lapply(group_files, function(file) {
    load(file)
    print(paste("Loaded:", file))
    
    # Loop through replicates and format
    replicates <- lapply(data, function(replicate) {
      replicate <- replicate %>% 
        select(t, Threshold)
      return(replicate)
    })
    replicates <- do.call('rbind', replicates)
    replicates <- replicates %>% 
      group_by(t) %>% 
      mutate(threshold = as.integer(as.character(cut(Threshold, 
                       breaks = c(-Inf, seq(0, 100, 0.25)), 
                       labels = seq(0, 100, 0.25))))) %>% 
      group_by(threshold) %>% 
      mutate(time = as.integer(as.character(cut(t, 
                     breaks = c(seq(0, 50000, 100), Inf),
                     labels = seq(0, 50000, 100), 
                     right = FALSE)))) %>% 
      group_by(time, threshold) %>% 
      summarise(count = length(threshold))
    return(replicates)
  })
  
  # Bind
  group_data <- do.call('rbind', group_data)
  all_values <- expand.grid(time = unique(group_data$time), threshold = unique(group_data$threshold))
  group_data <- group_data %>% 
    group_by(time, threshold) %>% 
    summarise(count = sum(count)) %>% 
    merge(all_values, by = c("time", "threshold"), all = TRUE)
  group_data$count[is.na(group_data$count)] <- 0
  group_data$freq <- group_data$count/max(group_data$count)
  
  # Plot
  myPalette <- colorRampPalette(brewer.pal(9, "PuBu"))
  # myPalette <- viridis_pal(option = "A", direction = -1)
  end_col <- viridis(length(group_sizes))[which(group_sizes == size)]
  max_color <- 0.03
  gg_threshtime <- ggplot(data = group_data,
                          aes(x = time, y = threshold, fill = freq, colour = freq)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradientn(colours = c("white", myPalette(100)),
                         limits = c(0, max_color),
                         oob = squish) +
    scale_colour_gradientn(colours =c("white", myPalette(100)),
                           limits = c(0, max_color),
                           oob = squish) +
    # scale_fill_gradientn(colours = c("white", end_col),
    #                      limits = c(0, max_color),
    #                      oob = squish) +
    # scale_colour_gradientn(colours = c("white", end_col),
    #                        limits = c(0, max_color),
    #                        oob = squish) +
    scale_x_continuous(name = "Time step",
                       breaks = seq(0, 50000, 10000),
                       label = comma,
                       expand = c(0, 0)) +
    scale_y_continuous(name = "Threshold",
                       breaks = seq(0, 100, 25),
                       limits = c(0, 100),
                       labels = c("0", "", "50", "", "100"),
                       expand = c(0, 0)) +
    theme(axis.text = element_text(colour = "black", size = 6),
          axis.title = element_text(size = 7),
          legend.position = "none",
          axis.ticks = element_line(size = 0.2, color = "black"),
          panel.border = element_rect(fill = NA, size = 0.3, color = "black"),
          panel.grid = element_blank(),
          plot.margin = unit(c(0.2, 0.35, 0.1, 0.1), "cm"))
  ggsave(gg_threshtime, file = paste0("output/ThresholdTime/ThreshTime_", size, ".png"), width = 90, height = 25, units = "mm", dpi = 600)
  ggsave(gg_threshtime, file = paste0("output/ThresholdTime/ThreshTime_", size, ".pdf"), width = 90, height = 25, units = "mm", dpi = 600)
  # smoothScatter(x = group_data$t, y = group_data$Threshold,
  #               ylim = c(0, 100))
  # plot(hexbin(x = group_data$t, y = group_data$Threshold, xbins = 30))
})







