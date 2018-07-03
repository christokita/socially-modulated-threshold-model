################################################################################
#
# Comparing various specialization plots
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


directory_path <- "output/Rdata/Sigma0-Epsilon0.1-Beta1.1/"
output_path <- "output/Rdata/_ProcessedData/"
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
list_folders <- folders[folders %in% c("Graphs", "Thresh", "Stim", "TaskTally", "Thresh1Time", "Thresh2Time")]

####################
# Load, bind, and save data
####################
for (folder in bind_folders) {
  files <- list.files(paste0(directory_path, folder), full.names = TRUE)
  for (file in files) {
    load(file)
    data <- do.call('rbind', data)
    if (!exists("compiled_data")) {
      compiled_data <- data
    } else {
      compiled_data <- rbind(compiled_data, data)
    }
  }
  compiled_data <- as.data.frame(compiled_data)
  save(compiled_data, file = paste0(output_path, folder, "/", run_info, ".Rdata"))
  rm(compiled_data)
}


####################
# Load, compile (lists), and save data
####################
for (folder in list_folders) {
  # Create subfolder to story data by group size
  full_output_path <- paste0(output_path, folder, "/", run_info, "/")
  dir.create(full_output_path, showWarnings = FALSE)
  # Get files
  files <- list.files(paste0(directory_path, folder), full.names = TRUE)
  # Get group sizes
  group_sizes <- gsub(".*/([0-9]+)-[0-9]+.Rdata", "\\1", files, perl = TRUE)
  group_sizes <- unique(group_sizes)
  # Loop through group sizes
  # Graphs are stored with different name
  if (folder == "Graphs") {
    for (i in group_sizes) {
      group_files <- files[grepl(paste0(".*/", i, "-[0-9]+.Rdata"), files)]
      # Loop through the files for this group size and bind
      for (file in group_files) {
        # Load data
        load(file)
        # Bind together
        ifelse(!exists("listed_data"),  listed_data <- ens_graphs, listed_data <- c(listed_data, ens_graphs))
      }
      save(listed_data, file = paste0(full_output_path, i, ".Rdata"))
      rm(listed_data)
    }
  # Check if threshold time series (massive so just move, don't process)
  } else if (folder %in% c("Thresh1Time", "Thresh2Time")) {
    file.copy(from = files, to = full_output_path)
    # Otherwise bind together normally  
  } else {
    for (i in group_sizes) {
      group_files <- files[grepl(paste0(".*/", i, "-[0-9]+.Rdata"), files)]
      # Loop through the files for this group size and bind
      for (file in group_files) {
        # Load data
        load(file)
        # Bind together
        ifelse(!exists("listed_data"),  listed_data <- data, listed_data <- c(listed_data, data))
      }
      save(listed_data, file = paste0(full_output_path, i, ".Rdata"))
      rm(listed_data)
    }
  }
}
