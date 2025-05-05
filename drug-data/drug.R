# ------------------------------------------------------------
# Clean environment and set working directory
# ------------------------------------------------------------
rm(list = ls())
setwd("C:/Users/User/Dropbox/PAPERS/projects/eleni/drug")

# ------------------------------------------------------------
# Load required libraries and helper functions
# ------------------------------------------------------------
library(igraph)
source("C:/Users/User/Dropbox/PAPERS/projects/eleni/drug/r_functions.R")
source("C:/Users/User/Dropbox/PAPERS/projects/eleni/drug/helper_functions_data_analysis.R")

# ------------------------------------------------------------
# Define number of nodes
# ------------------------------------------------------------
I <- 193

# ------------------------------------------------------------
# Load edge list and construct adjacency matrix
# ------------------------------------------------------------
dat <- as.matrix(read.csv("drug.txt", header = FALSE))
Ycube <- matrix(0, I, I)
for (k in 1:nrow(dat)) {
     i <- dat[k, 1]
     j <- dat[k, 2]
     Ycube[i, j] <- 1
     Ycube[j, i] <- 1
}

# ------------------------------------------------------------
# Basic structure checks
# ------------------------------------------------------------
isSymmetric(Ycube)        # Check if adjacency matrix is symmetric
any(rowSums(Ycube) == 0)  # Check for isolated nodes

# ------------------------------------------------------------
# Vectorize adjacency matrix and count total edges
# ------------------------------------------------------------
Y <- get_vec_adjacency(Ycube)
sum(Y)                    # Total number of edges (undirected)

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
pdf(file = "drug_sociomatrix.pdf")
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
plot_adjacency(Y = Ycube, labels = clust$membership)
dev.off()

# ------------------------------------------------------------
# Save graph plot
# ------------------------------------------------------------
pdf(file = "drug_graph.pdf")
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

# ------------------------------------------------------------
# Save processed data
# ------------------------------------------------------------
save(I, Y, Ycube, file = "drug.RData")