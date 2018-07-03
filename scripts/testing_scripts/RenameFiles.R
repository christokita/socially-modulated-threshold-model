
rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

run <- "Sigma0-Epsilon0.1-Beta1.1"

####################
# Go through folders and pad file names
####################
# Divide up folders into those where bound dataframes will be made
# and those where compiled lists will be made (for the sake of memory)
list_folders <- c("Graphs", "Thresh", "Stim", "TaskTally")

for (folder in list_folders) {
  print(folder)
  # Get directory name
  directory <- paste0("output/Rdata/_ProcessedData/", folder, "/", run, "/")
  # Get list of files
  files <- list.files(directory, full.names = TRUE)
  if (length(files > 0)) {
    # Create list of new file names
    new_names <- gsub(".*/([0-9]+)\\.Rdata", "\\1", files, perl = T)
    new_names <- str_pad(new_names, 3, pad = "0")
    new_names <- paste0(directory, new_names, ".Rdata")
    # Rename
    file.rename(from = files, to = new_names)
  }
}

other_folders <- c("Thresh1Time", "Thresh2Time")

for (folder in other_folders) {
  print(folder)
  # Get directory name
  directory <- paste0("output/Rdata/_ProcessedData/", folder, "/", run, "/")
  # Get list of files
  files <- list.files(directory, full.names = TRUE)
  if (length(files > 0)) {
    # Create list of new file names
    new_names <- gsub(".*/([0-9]+)-[0-9]+\\.Rdata", "\\1", files, perl = T)
    chunks <- gsub(".*/[0-9]+-([0-9]+)\\.Rdata", "\\1", files, perl = T)
    new_names <- str_pad(new_names, 3, pad = "0")
    new_names <- paste(new_names, chunks, sep = "-")
    new_names <- paste0(directory, new_names, ".Rdata")
    # Rename
    file.rename(from = files, to = new_names)
  }
}
