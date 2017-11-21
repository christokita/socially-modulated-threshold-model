################################################################################
#
# Examine Network from last run
#
################################################################################

g <- g_tot / gens

# Drop edges below baseline interaction probability
g[g < p] <- 0

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

write.csv(edgelist, file = "output/Networks/GroupSize70edgelistNoBias.csv", row.names = FALSE)
write.csv(nodelist, file = "output/Networks/GroupSize70nodelistNoBias.csv", row.names = FALSE)
