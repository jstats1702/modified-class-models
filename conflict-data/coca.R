# Load required packages
suppressMessages(suppressWarnings({
  library(igraph)
  library(ggplot2)
  library(corrplot)
}))

setwd("~/Dropbox/PAPERS/projects/eleni/conflict")

load("100he_2casos_organizaciones.RData")

source("~/Dropbox/PAPERS/projects/eleni/conflict/r_functions.R")
source("~/Dropbox/PAPERS/projects/eleni/conflict/helper_functions_data_analysis.R")

g <- metrics_list$graph_complete

# Remove names and isolate nodes
g <- delete_vertices(g, which(degree(g) == 0))

# Updated adjacency matrix without isolated nodes
Y <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))

rownames(Y) <- NULL
colnames(Y) <- NULL

g <- graph_from_adjacency_matrix(Y, mode = "undirected", diag = FALSE)

# Community detection
cl <- cluster_fast_greedy(g)
K <- length(unique(membership(cl)))

# Save data
Ycube <- Y
I     <- dim(Ycube)[1]
Y     <- get_vec_adjacency(Ycube)

save(I, Ycube, Y, file = "conflict.RData")

# ------------------------------------------------------------
# Community detection via fast greedy algorithm
# ------------------------------------------------------------
clust <- cluster_fast_greedy(g)
table(clust$membership)

# ------------------------------------------------------------
# Save sociomatrix plot
# ------------------------------------------------------------
pdf(file = "conflict_sociomatrix.pdf")
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
plot_adjacency(Y = Ycube, labels = clust$membership)
dev.off()

# ------------------------------------------------------------
# Save graph plot
# ------------------------------------------------------------
pdf(file = "conflict_graph.pdf")
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