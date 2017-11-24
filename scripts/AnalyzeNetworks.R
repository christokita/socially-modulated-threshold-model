################################################################################
#
# Analyze Networks of Social Interactions
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


filename <- "Sigma001-Eps001-Phi001-ConnectP01-Bias1.1"

####################
# Compare entropies
####################
# Load social
load("output/Rdata/Sigma001-Eps001-Phi001-ConnectP01-Bias1.1.Rdata")

number <- 180

g <- unlist(groups_graphs, recursive = FALSE)
g <- g[[number]]

degree <- rowSums(g)
degree <- as.data.frame(degree)
degree$Id <- row.names(degree)

threshMat <- unlist(groups_thresh, recursive = FALSE)
threshMat <- threshMat[[number]]

thresh <- threshMat %>% 
  as.data.frame(.) %>% 
  mutate(ThreshRatio = log(Thresh1 / Thresh2),
         Id = row.names(.))

X_tot <- unlist(groups_taskDist, recursive = FALSE)
X_tot <- X_tot[[number]]

activity <- X_tot %>% 
  as.data.frame(.) %>% 
  mutate(ActRatio = log(Task1 / Task2),
         ActTotal = Task1 + Task2,
         Id = row.names(.))

mergedNodes <- merge(degree, thresh)
mergedNodes <- merge(mergedNodes, activity)

qplot(data = mergedNodes, x = ThreshRatio, y = ActTotal, colour = degree) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, fill = NA))

