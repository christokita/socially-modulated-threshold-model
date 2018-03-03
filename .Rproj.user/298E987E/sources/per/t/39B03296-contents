################################################################################
#
# Examine Network from last run
#
################################################################################
rm(list = ls())
source("scripts/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


# Load social
load("output/Rdata/Sigma0.01-FIXED-Bias1.1.Rdata")


g <- groups_graphs[[6]][[1]]
threshMat <- groups_thresh[[6]][[1]]


diag(g) <- NA
# Only show edges above 50th percentile
percentiles <- quantile(g, na.rm = TRUE)
fiftypercent <- percentiles[3]
seventyfivepercent <- percentiles[4]
g[g < fiftypercent] <- 0
diag(g) <- 0

# Edgelist
graph <- graph.adjacency(g, weighted = T)
edgelist <- get.edgelist(graph)
edgelist <- as.data.frame(edgelist)
names(edgelist) <- c("Source", "Target")
edgelist$Weight <- E(graph)$weight 

# Nodelist
nodelist <- data.frame(Id = rownames(threshMat))
nodelist$Thresh1 <- threshMat[ , 1]
nodelist$Thresh2 <- threshMat[ , 2]
nodelist$ThreshRatio <- log(threshMat[ , 1] / threshMat[ , 2])
nodelist$ThreshRatioBounded <- nodelist$ThreshRatio
nodelist$ThreshRatioBounded[nodelist$ThreshRatioBounded < -0.5] <- -0.5 
nodelist$ThreshRatioBounded[nodelist$ThreshRatioBounded > 0.5] <- 0.5 

write.csv(edgelist, file = "output/Networks/GroupSize50FIXEDedgelist.csv", row.names = FALSE)
write.csv(nodelist, file = "output/Networks/GroupSize50FIXEDnodelist.csv", row.names = FALSE)

# Try plotting with igraph
library(igraph)
library(RColorBrewer)
g <- graph_from_adjacency_matrix(g, mode = "undirected", weighted = TRUE)
V(g)$ThreshRatio <- nodelist$ThreshRatio
E(g)$Reweight <- (E(g)$weight - min(E(g)$weight)) / (max(E(g)$weight) - min(E(g)$weight))
my.col <- colorRampPalette(brewer.pal(11, "RdYlBu"))

graphCol = palette(fine)[as.numeric(cut(btw,breaks = fine))]

plot(g, vertex.label = NA, edge.width = (E(g)$Reweight+0.01) * 5, vertex.color = V(g)$ThreshRatio)

# Try plotting with ggnet
library(GGally)

g_net <- network(g, directed = F)
ggnet2(g_net, color = )
