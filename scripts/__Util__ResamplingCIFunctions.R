##################################################
#
# Resampling CI Functions
#
##################################################


####################
# Generate CI
####################
resampleCI <- function(TaskDistList, resampleNumber) {
  # Get last group size (this will be resampled from)
  lastSize <- length(TaskDistList)
  toResample <- TaskDistList[[lastSize]]
  # Get group sizes of smaller groups
  sizes <- lapply(1:length(TaskDistList), function(x) {
    set <- TaskDistList[[x]]
    len <- nrow(set[[1]])
    return(len)
  })
  sizes <- unlist(sizes)
  # Resample
  resampledStats <- resampleStatistics(sizes = sizes, largestGroupList = toResample, resampleNumber = resampleNumber)
  resampledStats <- do.call("rbind", resampledStats)
  # Calculate mean value and sd by size
  sumStats <- resampledStats %>% 
    group_by(n) %>% 
    summarise(Mean1_avg = mean(Mean1),
              Mean1_95CI  = 1.96 * (sd(Mean1) / sqrt(length(Mean1))),
              SD1_avg = mean(SD1),
              SD1_95CI = 1.96 * (sd(SD1) / sqrt(length(SD1))),
              SD2_avg = mean(SD2),
              SD2_95CI = 1.96 * (sd(SD2) / sqrt(length(SD2))),
              SD_avg = mean(SD),
              SD_95CI = 1.96 * (sd(SD) / sqrt(length(SD))),
              SD_95Dist = 1.96 * sd(SD))
  return(sumStats)
}


####################
# Resample values from the largest group size
####################
resampleStatistics <- function(sizes, largestGroupList, resampleNumber) {
  resampledStats <- lapply(sizes, function(size) {
    # Repeat the resampling process
    replicateSamp <- list()
    for (i in 1:resampleNumber) {
      # Loop through replicates in largest group size
      sizeStats <- lapply(1:length(largestGroupList), function(x) {
        rep <- largestGroupList[[x]]
        # Sample required number from each colony
        sampledDist1 <- sample(rep$Task1, size = size, replace = FALSE)
        sampleSD1 <- sd(sampledDist1)
        sampleMean1 <- mean(sampledDist1)
        sampledDist2 <- sample(rep$Task2, size = size, replace = FALSE)
        sampleSD2 <- sd(sampledDist2)
        sampleMean2 <- mean(sampledDist2)
        # Put in dataframe and return
        sampledDist <- data.frame(SD1 = sampleSD1, Mean1 = sampleMean1,
                                  SD2 = sampleSD2, Mean2 = sampleMean2)
        return(sampledDist)
      })
      # Unlist and bind
      sizeStats <- do.call("rbind", sizeStats)
      # Get mean value of SD and Mean
      stats <- sizeStats %>% 
        summarise(n = size, SD1 = mean(SD1), Mean1 = mean(Mean1), SD2 = mean(SD2), Mean2 = mean(Mean2)) %>% 
        mutate(SD = (SD1 + SD2) / 2)
      replicateSamp[[i]] <- stats
    }
    # Bind
    replicateSamp <- do.call("rbind", replicateSamp)
    return(replicateSamp)
  })
}


####################
# Generate CI
####################
resampleCI_experiment <- function(ExperimentalData, resampleNumber) {
  # Get last group size (this will be resampled from)
  lastSize <- max(ExperimentalData$GS)
  toResample <- ExperimentalData %>% 
    filter(GS == lastSize)
  # Get group sizes of smaller groups
  sizes <- unique(ExperimentalData$GS)
  # Resample
  resampledStats <- resampleStatistics_experiment(sizes = sizes, largestGroupList = toResample, resampleNumber = resampleNumber)
  resampledStats <- do.call("rbind", resampledStats)
  # Calculate mean value and sd by size
  sumStats <- resampledStats %>% 
    group_by(n) %>% 
    summarise(Mean_avg = mean(Mean),
              Mean_95CI  = 1.96 * (sd(Mean) / sqrt(length(Mean))),
              Mean_95Dist = sd(Mean) * 1.96,
              SD_avg = mean(SD),
              SD_95Dist = sd(SD) * 1.96,
              SD_95CI = 1.96 * (sd(SD) / sqrt(length(SD))))
  return(sumStats)
}


####################
# Resample values from the largest group size
####################
resampleStatistics_experiment <- function(sizes, largestGroupList, resampleNumber) {
  # Unique colonies
  largestGroupList <- largestGroupList %>% 
    filter(!is.na(RMSDnorm))
  colonies <- unique(largestGroupList$boxID)
  # sample
  resampledStats <- lapply(sizes, function(size) {
    # Repeat the resampling process
    replicateSamp <- list()
    for (i in 1:resampleNumber) {
      # Loop through replicates in largest group size
      sizeStats <- lapply(1:length(colonies), function(x) {
        rep <- largestGroupList[largestGroupList$boxID == colonies[x], ]
        # Sample required number from each colony
        sampledDist <- sample(rep$RMSDnorm, size = size, replace = FALSE)
        sampleSD <- sd(sampledDist)
        sampleMean <- mean(sampledDist)
        # Put in dataframe and return
        sampledDist <- data.frame(SD = sampleSD, Mean = sampleMean)
        return(sampledDist)
      })
      # Unlist and bind
      sizeStats <- do.call("rbind", sizeStats)
      # Get mean value of SD and Mean
      stats <- sizeStats %>% 
        summarise(n = size, SD = mean(SD), Mean = mean(Mean))
      replicateSamp[[i]] <- stats
    }
    # Bind
    replicateSamp <- do.call("rbind", replicateSamp)
    return(replicateSamp)
  })
}

