################################################################################
#
# Analyze stimulus data
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")

####################
# Average change in stim levels: beginning to end
####################
# Social thresholds
files <- list.files("output/Rdata/_ProcessedData/Stim/Sigma0-Epsilon0.1-Beta1.1/", full.names = TRUE)
stim_data <- lapply(files, function(file) {
  # Load group size data
  load(file)
  n <- as.numeric(gsub(".*/([0-9]+)\\.Rdata", "\\1", file, perl = T))
  print(paste("Loaded: Group Size", n))
  # Summarise within each replicate
  rep_data <- lapply(listed_data, function(replicate) {
    replicate <- as.data.frame(replicate)
    replicate <- replicate[-1, ] #remove time step 0
    replicate$sTotal <- replicate$s1 + replicate$s2
    # Grab first and last 1000 time steps
    begin <- head(replicate, 1000)
    end <- tail(replicate, 1000)
    # Summarise
    begin <- begin %>% 
      summarise(s1 = mean(s1),
                s2 = mean(s2),
                sTotal = mean(sTotal))
    end <- end %>% 
      summarise(s1 = mean(s1),
                s2 = mean(s2),
                sTotal = mean(sTotal))
    diff <- end - begin
    return(diff)
  })
  rep_data <- do.call('rbind', rep_data)
  # Combine all stim data
  all_stim <- c(rep_data$s1, rep_data$s2)
  # Calculate statistics and return
  to_return <- data.frame(n = n, 
                          sMean = mean(rep_data$sTotal), 
                          sSD = sd(rep_data$sTotal), 
                          sSE = sd(rep_data$sTotal)/sqrt(nrow(rep_data)))
  return(to_return)
})
# Bind
stim_data <- do.call("rbind", stim_data)

# Fixed thresholds
files <- list.files("output/Rdata/_ProcessedData/Stim/Sigma0-Epsilon0.1-Beta1.1/", full.names = TRUE)
stim_data <- lapply(files, function(file) {
  # Load group size data
  load(file)
  n <- as.numeric(gsub(".*/([0-9]+)\\.Rdata", "\\1", file, perl = T))
  print(paste("Loaded: Group Size", n))
  # Summarise within each replicate
  rep_data <- lapply(listed_data, function(replicate) {
    replicate <- as.data.frame(replicate)
    replicate <- replicate[-1, ] #remove time step 0
    replicate$sTotal <- replicate$s1 + replicate$s2
    # Grab first and last 1000 time steps
    begin <- head(replicate, 1000)
    end <- tail(replicate, 1000)
    # Summarise
    begin <- begin %>% 
      summarise(s1 = mean(s1),
                s2 = mean(s2),
                sTotal = mean(sTotal))
    end <- end %>% 
      summarise(s1 = mean(s1),
                s2 = mean(s2),
                sTotal = mean(sTotal))
    diff <- end - begin
    return(diff)
  })
  rep_data <- do.call('rbind', rep_data)
  # Combine all stim data
  all_stim <- c(rep_data$s1, rep_data$s2)
  # Calculate statistics and return
  to_return <- data.frame(n = n, 
                          sMean = mean(rep_data$sTotal), 
                          sSD = sd(rep_data$sTotal), 
                          sSE = sd(rep_data$sTotal)/sqrt(nrow(rep_data)))
  return(to_return)
})
# Bind
stim_data <- do.call("rbind", stim_data)

# Plot
gg_stimdiff <- ggplot(stim_data, aes(x = n, y = sMean)) +
  geom_hline(yintercept = 0, 
             color = "grey90", 
             size = 0.3) +
  geom_errorbar(aes(ymin = sMean - sSD, ymax = sMean + sSD), 
                size = 0.3,
                width = 0, 
                color = "#4d4d4d") +
  geom_point(size = 0.8, 
             color = "#4d4d4d", 
             fill = "#4d4d4d") +
  theme_classic() +
  ylab("Change in Total Stimulus") +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.line = element_line(size = 0.3, color = "black"),
        aspect.ratio = 1)

gg_stimdiff

ggsave("output/StimLevels/ChangeInStimLevels.png", width = 45, height = 45, units = 'mm', dpi = 400)


####################
# Average change in stim levels: beginning to end
####################
# Runs of intertest
runs <- c("Sigma0.05-Epsilon0-Beta1.1", 
          "Sigma0-Epsilon0.1-Beta1.1")
run_names <- c("Fixed", "Social")


stim_data_comp <- lapply(1:length(runs), function(run) {
  # Load stim data
  print(runs[run])
  model <- run_names[run]
  files <- list.files(paste0("output/Rdata/_ProcessedData/Stim/", runs[run], "/"), full.names = TRUE)
  run_data <- lapply(files, function(file) {
    load(file)
    n <- as.numeric(gsub(".*/([0-9]+)\\.Rdata", "\\1", file, perl = T))
    print(paste("Loaded: Group Size", n))
    # Summarise within each replicate
    rep_data <- lapply(listed_data, function(replicate) {
      replicate <- as.data.frame(replicate)
      replicate <- replicate[-1, ] #remove time step 0
      replicate$sTotal <- replicate$s1 + replicate$s2
      # Grab first and last 1000 time steps
      begin <- head(replicate, 1000)
      end <- tail(replicate, 1000)
      # Summarise
      begin <- begin %>% 
        summarise(s1 = mean(s1),
                  s2 = mean(s2),
                  sTotal = mean(sTotal))
      end <- end %>% 
        summarise(s1 = mean(s1),
                  s2 = mean(s2),
                  sTotal = mean(sTotal))
      diff <- end - begin
      return(diff)
    })
    # Bind, label, return by group size within run
    rep_data <- do.call("rbind", rep_data)
    rep_data$n <- n
    return(rep_data)
  })
  # Bind, label, return by run
  run_data <- do.call("rbind", run_data)
  run_data$Model <- model
  return(run_data)
})

# Bind and summarise
stim_data_comp <- do.call('rbind', stim_data_comp) %>% 
  group_by(Model, n) %>% 
  summarise(sMean = mean(sTotal), 
            sSD = sd(sTotal), 
            sSE = sd(sTotal)/sqrt(length(n)))

# Plot
gg_stimdiff_comp <- ggplot(stim_data_comp, aes(x = n, y = sMean, group = Model, fill = Model, color = Model)) +
  geom_hline(yintercept = 0, 
             color = "grey90", 
             size = 0.3) +
  geom_errorbar(aes(ymin = sMean - sSD, ymax = sMean + sSD), 
                size = 0.3,
                width = 0) +
  geom_point(size = 0.8,
             shape = 21) +
  theme_classic() +
  ylab("Change in Total Stimulus") +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_color_manual(name = "Threshold type",
                     values = c("#878787", "#4d4d4d")) +
  scale_fill_manual(name = "Threshold type",
                    values = c("#ffffff", "#4d4d4d")) +
  theme_ctokita() +
  theme(aspect.ratio = 1)
gg_stimdiff_comp

ggsave("output/StimLevels/ChangeInStimLevels_Comparison.png", height = 45, width = 75, units = 'mm', dpi = 800)

####################
# Look at individual examples
####################
# load(files[7])
# 
# example <- as.data.frame(listed_data[[1]])
# example <- example %>% 
#   mutate(TotalStim = s1 + s2)
# 
# gg_example <- ggplot(example, aes(x = t, y = TotalStim)) +
#   geom_line(size = 0.5) +
#   theme_classic() +
#   ylab("Change in Total Stim. Level") +
#   scale_x_continuous(breaks = seq(0, 50000, 10000)) +
#   theme(axis.text = element_text(colour = "black", size = 6),
#         axis.title = element_text(size = 7),
#         axis.ticks = element_line(size = 0.3, color = "black"),
#         axis.line = element_line(size = 0.3, color = "black"),
#         aspect.ratio = 1)
# 
# gg_example

####################
# Average time series by group size
####################
# files <- list.files("output/Rdata/_ProcessedData/Stim/Sigma0-Epsilon0.1-Beta1.1/", full.names = TRUE)
# stim_time <- lapply(files, function(file) {
#   # Load group size data
#   load(file)
#   n <- as.numeric(gsub(".*/([0-9]+)\\.Rdata", "\\1", file, perl = T))
#   print(paste("Loaded: Group Size", n))
#   # Summarise within each replicate
#   rep_data <- lapply(listed_data, function(replicate) {
#     replicate <- as.data.frame(replicate)
#     replicate <- replicate %>% 
#       mutate(sTotal = s1 + s2) %>% 
#       select(sTotal)
#     return(replicate)
#   })
#   # Bind and summarise
#   rep_data <- do.call('cbind', rep_data)
#   colnames(rep_data) <- paste0('sTotal', 1:ncol(rep_data))
#   to_return <- rep_data %>% 
#     mutate(sTotal = rowSums(.) / 100) %>% 
#     select(sTotal) %>% 
#     mutate(t = 0:(nrow(.)-1)) %>% 
#     mutate(n = n) %>% 
#     select(n, t, sTotal)
#   # Combine all stim data
#   return(to_return)
# })
# # Bind
# stim_time <- do.call("rbind", stim_time)
# 
# # Plot
# select_sizes <- stim_time %>% 
#   filter(n %in% seq(70, 100)) %>%
#   mutate(n = as.factor(n))
# 
# pal <- brewer.pal(9, "PuRd")
# 
# gg_stimtime <- ggplot(select_sizes, aes(x = t, y = sTotal, group = n, color = n)) +
#   geom_line(size = 0.5) +
#   theme_classic() +
#   ylab("Total stimulus") +
#   scale_x_continuous(breaks = seq(0, 50000, 10000)) +
#   scale_color_manual(values = pal) +
#   theme(axis.text = element_text(colour = "black", size = 6),
#         axis.title = element_text(size = 7),
#         axis.ticks = element_line(size = 0.3, color = "black"),
#         axis.line = element_line(size = 0.3, color = "black"),
#         aspect.ratio = 1)
# 
# gg_stimtime
# 
# ggsave("output/StimLevels/TimeSeries.png", width = 180, height = 180, units = "mm", dpi = 800)
