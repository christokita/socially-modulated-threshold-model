################################################################################
#
# Ensure that rank corr run and original run are indeed equal so that
# Rank corr data can simply be added to the original run
#
################################################################################

##########################
# Comparing threshold distirbutions
##########################

##### Rank Corr Run #####
files <- list.files(paste0("output/Rdata/archive/Thresh/Sigma0-Epsilon0.1-Beta1.1_RankCorr/"), full.names = TRUE)
# Loop through group sizes
group_thresh <- lapply(files, function(file) {
  # Load
  load(file)
  # Bind
  thresh_data <- as.data.frame(do.call("rbind", listed_data))
  return(thresh_data)
})

# Bind 
rank_data <- do.call('rbind', group_thresh)

##### Original Run #####
files <- list.files(paste0("output/Rdata/_ProcessedData/Thresh/Sigma0-Epsilon0.1-Beta1.1/"), full.names = TRUE)
# Loop through group sizes
group_thresh <- lapply(files, function(file) {
  # Load
  load(file)
  # Bind
  thresh_data <- as.data.frame(do.call("rbind", listed_data))
  return(thresh_data)
})

# Bind 
original_data <- do.call('rbind', group_thresh)

# Test
table(original_data == rank_data)


##########################
# Comparing entropy data
##########################

rm(list = ls())

##### Rank Corr Run #####
load("output/Rdata/archive/Entropy/Sigma0-Epsilon0.1-Beta1.1_RankCorr.Rdata")
rank_data <- compiled_data %>% 
  arrange(n, sim, chunk)

##### Original Run #####
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")
original_data <- compiled_data  %>% 
  arrange(n, sim, chunk)

# Test
table(original_data == rank_data)
