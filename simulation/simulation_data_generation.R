# ------------------------------------------------------------
# Setup: working directory and clean environment
# ------------------------------------------------------------
setwd("C:/Users/User/Dropbox/PAPERS/projects/eleni")
rm(list = ls())

# ------------------------------------------------------------
# Load required libraries and custom functions
# ------------------------------------------------------------
library(igraph)
library(corrplot)

source("C:/Users/User/Dropbox/PAPERS/projects/eleni/r_functions.R")
source("C:/Users/User/Dropbox/PAPERS/projects/eleni/helper_functions_simulation.R")

# ------------------------------------------------------------
# Define simulation parameters
# ------------------------------------------------------------
p_within_values  <- seq(0.50, 0.60, by = 0.02)
p_between_values <- seq(0.10, 0.15, by = 0.01)
cluster_type <- "random"

# ------------------------------------------------------------
# Simulate SBM networks with fixed seed
# ------------------------------------------------------------
set.seed(42)
sbm_simulation_1 <- generate_SBM_network(
     num_nodes = 100, num_communities = 5,
     p_within_values, p_between_values, 
     cluster_type
)

set.seed(42)
sbm_simulation_2 <- generate_SBM_network(
     num_nodes = 200, num_communities = 10,
     p_within_values, p_between_values, 
     cluster_type
)

# ------------------------------------------------------------
# Reorder adjacency matrices for plotting
# ------------------------------------------------------------
tmp_1 <- ordering0(AM = sbm_simulation_1$adjacency_matrix, labels = sbm_simulation_1$communities)
tmp_2 <- ordering0(AM = sbm_simulation_2$adjacency_matrix, labels = sbm_simulation_2$communities)

# ------------------------------------------------------------
# Save sociomatrix plot for scenario 1
# ------------------------------------------------------------
pdf(file = "scenario_1_sociomatrix.pdf")
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
plot_adjacency(Y = tmp_1$AM, labels = tmp_1$labels)
dev.off()

# ------------------------------------------------------------
# Save graph plot for scenario 1
# ------------------------------------------------------------
pdf(file = "scenario_1_graph.pdf")
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
set.seed(123)
plot(sbm_simulation_1$graph,
     layout = layout_with_fr,
     mark.groups = NULL,
     edge.color = adjustcolor("black", 0.15),
     vertex.size = 3,
     vertex.color = sbm_simulation_1$communities,
     vertex.frame.color = sbm_simulation_1$communities,
     vertex.label = NA)
dev.off()

# ------------------------------------------------------------
# Save sociomatrix plot for scenario 2
# ------------------------------------------------------------
pdf(file = "scenario_3_sociomatrix.pdf")
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
plot_adjacency(Y = tmp_2$AM, labels = tmp_2$labels)
dev.off()

# ------------------------------------------------------------
# Save graph plot for scenario 2
# ------------------------------------------------------------
pdf(file = "scenario_2_graph.pdf")
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
set.seed(123)
plot(sbm_simulation_2$graph,
     layout = layout_with_fr,
     mark.groups = NULL,
     edge.color = adjustcolor("black", 0.15),
     vertex.size = 3,
     vertex.color = sbm_simulation_2$communities,
     vertex.frame.color = sbm_simulation_2$communities,
     vertex.label = NA)
dev.off()
