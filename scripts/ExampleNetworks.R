################################################################################
#
# Examine Network from last run
#
################################################################################

g <- g_tot / gens
diag(g) <- NA
# Only show edges above 50th percentile
percentiles <- quantile(g, na.rm = TRUE)
fiftypercent <- percentiles[3]
g[g <= fiftypercent] <- 0
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

write.csv(edgelist, file = "output/Networks/GroupSize20edgelist.csv", row.names = FALSE)
write.csv(nodelist, file = "output/Networks/GroupSize20nodelist.csv", row.names = FALSE)

# Try plotting with igraph
library(igraph)
library(RColorBrewer)
g <- graph_from_adjacency_matrix(g, mode = "undirected", weighted = TRUE)
V(g)$ThreshRatio <- nodelist$ThreshRatio
E(g)$Reweight <- (E(g)$weight - min(E(g)$weight)) / (max(E(g)$weight) - min(E(g)$weight))
my.col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(diff(range(V(g)$ThreshRatio)))
plot(g, vertex.label = NA, edge.width = (E(g)$Reweight+0.01) * 5, vertex.color = V(g)$ThreshRatio)
