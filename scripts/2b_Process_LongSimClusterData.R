################################################################################
#
# Comparing various specialization plots
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")


directory_path <- "output/Rdata/long_sims/LongSim-Sigma0-Epsilon0.5-Beta1.05/"
output_path <- "output/Rdata/_ProcessedData/_LongSims/"
run_info <- gsub("^.*(Sigma.*)/$", "\\1", directory_path, perl = TRUE)

####################
# Create folders
####################
# Get names of folders
folders <- list.files(directory_path)
output_folders <- list.files(output_path)

# Create folders for those not existing in processed data
missing_folders <- folders[!folders %in% output_folders]
for (missing_folder in missing_folders) {
  dir.create(paste0(output_path, missing_folder))
}

# Divide up folders into those where bound dataframes will be made
# and those where compiled lists will be made (for the sake of memory)
bind_folders <- folders[folders %in% c("Entropy", "TaskDist")]
list_folders <- folders[folders %in% c("Graphs", "Thresh", "Stim")]

####################
# Bind and save entropy
####################
files <- list.files(paste0(directory_path, "Entropy"), full.names = T)
entropy_data <- list()
for (i in 1:length(files)) {
  file <- files[i]
  sim_number <- as.numeric(gsub(file, pattern = ".*Sim_([0-9]+).Rdata", replacement = "\\1", perl = T))
  load(file)
  entropy <- as.data.frame(entropy)
  entropy$sim <- sim_number
  entropy_data[[i]] <- entropy
} 
entropy_data <- do.call("rbind", entropy_data)
save(entropy_data, file = paste0(output_path, "/Entropy/", run_info, ".Rdata"))

####################
# Bind and save task distributions
####################
files <- list.files(paste0(directory_path, "TaskDist"), full.names = T)
taskdist_data <- list()
for (i in 1:length(files)) {
  file <- files[i]
  sim_number <- as.numeric(gsub(file, pattern = ".*Sim_([0-9]+).Rdata", replacement = "\\1", perl = T))
  load(file)
  totalTaskDist <- as.data.frame(totalTaskDist)
  totalTaskDist$sim <- sim_number
  taskdist_data[[i]] <- totalTaskDist
} 
taskdist_data <- do.call("rbind", taskdist_data)
save(taskdist_data, file = paste0(output_path, "/TaskDist/", run_info, ".Rdata"))

####################
# Bind and save task distributions
####################
files <- list.files(paste0(directory_path, "TaskDist"), full.names = T)
taskdist_data <- list()
for (i in 1:length(files)) {
  file <- files[i]
  sim_number <- as.numeric(gsub(file, pattern = ".*Sim_([0-9]+).Rdata", replacement = "\\1", perl = T))
  load(file)
  totalTaskDist <- as.data.frame(totalTaskDist)
  totalTaskDist$sim <- sim_number
  taskdist_data[[i]] <- totalTaskDist
} 
taskdist_data <- do.call("rbind", taskdist_data)
save(taskdist_data, file = paste0(output_path, "/TaskDist/", run_info, ".Rdata"))

####################
# Bind and save graphs
####################
files <- list.files(paste0(directory_path, "Graphs"), full.names = T)
graphs_data <- list()
for (i in 1:length(files)) {
  file <- files[i]
  sim_number <- as.numeric(gsub(file, pattern = ".*Sim_([0-9]+).Rdata", replacement = "\\1", perl = T))
  load(file)
  graphs_data[[i]] <- g_tot
} 
save(graphs_data, file = paste0(output_path, "/Graphs/", run_info, ".Rdata"))


