# ------------------------------------------------------------
# Clean environment and set working directory
# ------------------------------------------------------------
rm(list = ls())
setwd("C:/Users/User/Dropbox/PAPERS/projects/eleni/salter")

# ------------------------------------------------------------
# Load required libraries and helper functions
# ------------------------------------------------------------
library(igraph)
source("C:/Users/User/Dropbox/PAPERS/projects/eleni/salter/r_functions.R")
source("C:/Users/User/Dropbox/PAPERS/projects/eleni/salter/helper_functions_data_analysis.R")

# ------------------------------------------------------------
# Load edge list and construct adjacency matrix
# ------------------------------------------------------------
load("C:/Users/User/Dropbox/PAPERS/projects/eleni/salter/salter.RData")

isolated_nodes <- colSums(Ycube) == 0
Ycube <- Ycube[!isolated_nodes, !isolated_nodes]

# ------------------------------------------------------------
# Basic structure checks
# ------------------------------------------------------------
isSymmetric(Ycube)        # Check if adjacency matrix is symmetric
any(rowSums(Ycube) == 0)  # Check for isolated nodes

# ------------------------------------------------------------
# Build igraph object and visualize raw graph
# ------------------------------------------------------------
g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")

# ------------------------------------------------------------
# Community detection via fast greedy algorithm
# ------------------------------------------------------------
clust <- cluster_fast_greedy(g)
table(clust$membership)

# ------------------------------------------------------------
# Save sociomatrix plot
# ------------------------------------------------------------
pdf(file = "salter_sociomatrix.pdf")
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
plot_adjacency(Y = Ycube, labels = clust$membership)
dev.off()

# ------------------------------------------------------------
# Save graph plot
# ------------------------------------------------------------
pdf(file = "salter_graph.pdf")
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
set.seed(123)
plot(g,
     layout = layout_with_fr,
     mark.groups = NULL,
     edge.color = adjustcolor("black", 0.05),
     vertex.size = 4,
     vertex.color = clust$membership,
     vertex.frame.color = clust$membership,
     vertex.label = NA)
dev.off()