class_dm_test <- function(Ycube, I, K, samples) {
     # Load necessary libraries and compiled C++ functions
     Rcpp::sourceCpp("class_dm_functions.cpp")
     
     # Compute observed statistics from the adjacency matrix
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(
          edge_density(g),
          transitivity(g, type = "global"),
          assortativity_degree(g),
          mean_distance(g, directed = FALSE),
          mean(degree(g)),
          sd(degree(g))
     )
     
     # Preallocate matrix for test statistics
     B <- nrow(samples$Xi_chain)
     n_test_stats <- length(test_obs)
     test_stats <- matrix(NA, B, n_test_stats)
     
     # Simulation loop
     for (b in seq_len(B)) {
          # Extract parameters
          Lambda <- samples$Lambda_chain[b, ]
          Xi     <- samples$Xi_chain[b, ]
          
          # Simulate adjacency matrix and convert to graph
          A <- simulate_data(I, K, Lambda, Xi)
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree and test statistics
          degrees <- degree(g)
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = FALSE),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat(sprintf("Simulation progress: %.1f%% completed\n", 100 * b / B))
          }
     }
     
     # Compute posterior predictive p-values (PPPs)
     ppps <- colMeans(test_stats < test_obs)
     
     # Return results as a named list
     list(
          test_obs   = test_obs,
          test_stats = test_stats,
          ppps       = ppps
     )
}

class_dp_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("class_dp_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Lambda <- samples$Lambda_chain[[b]]
          Xi     <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          A <- simulate_data(I, Lambda, Xi)
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

class_py_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("class_py_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Lambda <- samples$Lambda_chain[[b]]
          Xi     <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          A <- simulate_data(I, Lambda, Xi)
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

class_gn_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("class_gn_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Lambda <- samples$Lambda_chain[[b]]
          Xi     <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          A <- simulate_data(I, Lambda, Xi)
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

clist_dm_test <- function(Ycube, I, K, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("clist_dm_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Eta <- samples$Eta_chain[b,]
          Xi  <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          A <- simulate_data(I, K, Eta, Xi)
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

clist_dp_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("clist_dp_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Eta <- samples$Eta_chain[[b]]
          Xi  <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          repeat { 
               A <- simulate_data(I, Eta, Xi)
               if(!any(is.na(A))) {
                    break
               }
          }
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

clist_py_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("clist_py_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Eta <- samples$Eta_chain[[b]]
          Xi  <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          repeat { 
               A <- simulate_data(I, Eta, Xi)
               if(!any(is.na(A))) {
                    break
               }
          }
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

clist_gn_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("clist_gn_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Eta <- samples$Eta_chain[[b]]
          Xi  <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          repeat { 
               A <- simulate_data(I, Eta, Xi)
               if(!any(is.na(A))) {
                    break
               }
          }
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

cleigen_dm_test <- function(Ycube, I, K, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("cleigen_dm_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Eta <- samples$Eta_chain[b,]
          Xi  <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          A <- simulate_data(I, K, Eta, Xi)
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

cleigen_dp_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("cleigen_dp_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Eta <- samples$Eta_chain[[b]]
          Xi  <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          repeat { 
               A <- simulate_data(I, Eta, Xi)
               if(!any(is.na(A))) {
                    break
               }
          }
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

cleigen_py_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("cleigen_py_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Eta <- samples$Eta_chain[[b]]
          Xi  <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          repeat { 
               A <- simulate_data(I, Eta, Xi)
               if(!any(is.na(A))) {
                    break
               }
          }
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}

cleigen_gn_test <- function(Ycube, I, samples) 
{
     # settings
     require(igraph)
     Rcpp::sourceCpp("cleigen_gn_functions.cpp")
     
     # observed statistics
     g <- graph_from_adjacency_matrix(Ycube, mode = "undirected")
     test_obs <- c(edge_density(g),
                   transitivity(g, type = "global"),
                   assortativity_degree(g),
                   mean_distance(g, directed = F),
                   mean(degree(g)),
                   sd(degree(g)))
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 6
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          Eta <- samples$Eta_chain[[b]]
          Xi  <- samples$Xi_chain[b,]
          
          # Simulate adjacency matrix and convert to graph
          repeat { 
               A <- simulate_data(I, Eta, Xi)
               if(!any(is.na(A))) {
                    break
               }
          }
          g <- graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degree
          degrees <- degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               edge_density(g),
               transitivity(g, type = "global"),
               assortativity_degree(g),
               mean_distance(g, directed = F),
               mean(degrees),
               sd(degrees)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     ppps <- c(mean(test_stats[,1] < test_obs[1]),
               mean(test_stats[,2] < test_obs[2]),
               mean(test_stats[,3] < test_obs[3]),
               mean(test_stats[,4] < test_obs[4]),
               mean(test_stats[,5] < test_obs[5]),
               mean(test_stats[,6] < test_obs[6]))
     
     # return
     list(test_obs = test_obs,
          test_stats = test_stats,
          ppps = ppps)
}