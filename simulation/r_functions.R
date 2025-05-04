# Format a number to `k` decimal places
dec <- function(x, k) formatC(round(x, k), format = "f", digits = k)

# Procrustes alignment: Align Z to Z0
procus <- function(Z, Z0, K) {
     # Center Z by aligning its columns with Z0
     Z <- sweep(Z, 2, colMeans(Z)) + colMeans(Z0)
     
     # Compute alignment matrix
     A <- t(Z) %*% (Z0 %*% t(Z0)) %*% Z
     eA <- eigen(A, symmetric = TRUE)
     Ahalf <- eA$vectors[, 1:K] %*% diag(sqrt(eA$values[1:K])) %*% t(eA$vectors[, 1:K])
     
     # Align Z to Z0 using Procrustes transformation
     t(t(Z0) %*% Z %*% solve(Ahalf) %*% t(Z))
}

# Infer the number of nodes (I) from a binary adjacency matrix (Y)
get_I <- function(Y) (1 + sqrt(1 + 8 * nrow(Y))) / 2

# Convert (i, ii) pair to a flattened index for a matrix of size I
get_k <- function(i, ii, I) I * (I - 1) / 2 - (I - i + 1) * (I - i) / 2 + ii - i

# Convert (k, kk) pair to a flattened index for a diagonal matrix of size K
get_k_diag <- function(k, kk, K) K * (k - 1) + kk - (k - 1) * k / 2

# Generate fold assignments for cross-validation
get_folds <- function(M, J, L) {
     folds <- matrix(NA_integer_, nrow = M, ncol = J)
     
     for (j in seq_len(J)) {
          fold_assignments <- rep(seq_len(L), length.out = M)  # Ensures balanced assignment
          folds[, j] <- sample(fold_assignments, M, replace = FALSE)  # Shuffle assignments
     }
     
     folds
}

generate_SBM_network <- function(num_nodes, num_communities, p_within_values, p_between_values, cluster_type = "equal") 
{
     suppressMessages(suppressWarnings(library(igraph)))
     
     # Validate cluster type
     if (!cluster_type %in% c("equal", "random")) {
          stop("Invalid cluster_type. Use 'equal' or 'random'.")
     }
     
     # --- Assign Community Sizes ---
     if (cluster_type == "equal") {
          # Divide nodes equally and distribute remainder
          base_size <- num_nodes %/% num_communities
          remainder <- num_nodes %% num_communities
          
          community_sizes <- rep(base_size, num_communities)
          if (remainder > 0) {
               community_sizes[1:remainder] <- community_sizes[1:remainder] + 1
          }
          
     } else {  # Random-sized communities
          # Ensure each community gets at least one node
          community_sizes <- rep(1, num_communities)
          remaining_nodes <- num_nodes - num_communities
          
          if (remaining_nodes > 0) {
               raw_sizes <- runif(num_communities)  # Generate weights
               extra_sizes <- floor((raw_sizes / sum(raw_sizes)) * remaining_nodes)
               
               # Adjust if rounding error occurs
               while (sum(extra_sizes) != remaining_nodes) {
                    adjust_index <- sample(num_communities, 1)
                    diff <- remaining_nodes - sum(extra_sizes)
                    
                    if (diff > 0) {
                         extra_sizes[adjust_index] <- extra_sizes[adjust_index] + 1
                    } else if (diff < 0 && extra_sizes[adjust_index] > 0) {
                         extra_sizes[adjust_index] <- extra_sizes[adjust_index] - 1
                    }
               }
               
               community_sizes <- community_sizes + extra_sizes
          }
     }
     
     # Final validation to correct any rounding excess
     while (sum(community_sizes) > num_nodes) {
          community_index <- sample(num_communities, 1, T, community_sizes)
          community_sizes[community_index] <- community_sizes[community_index] - 1
     }
     
     # Assign communities
     node_communities <- rep(seq_len(num_communities), times = community_sizes)
     
     # Ensure correct number of nodes
     if (length(node_communities) != num_nodes) {
          stop("Error: node_communities does not match num_nodes!")
     }
     
     # --- Construct Probability Matrix ---
     P <- matrix(0, nrow = num_communities, ncol = num_communities)
     
     # Set within-community connection probabilities
     P[diag(num_communities) == 1] <- ifelse(length(p_within_values) == 1, 
                                             p_within_values, 
                                             sample(p_within_values, num_communities, replace = TRUE))
     
     # Set between-community probabilities (ensure symmetry)
     for (k in 1:(num_communities - 1)) {
          for (l in (k + 1):num_communities) {
               p_value <- ifelse(length(p_between_values) == 1, 
                                 p_between_values, 
                                 sample(p_between_values, 1))
               P[k, l] <- p_value
               P[l, k] <- p_value
          }
     }
     
     # --- Generate Adjacency Matrix ---
     adjacency_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
     
     for (i in 1:(num_nodes - 1)) {
          for (j in (i + 1):num_nodes) {
               k_i <- node_communities[i]
               k_j <- node_communities[j]
               p_edge <- P[k_i, k_j]
               
               # Sample edge presence
               edge <- rbinom(1, 1, p_edge)
               adjacency_matrix[i, j] <- edge
               adjacency_matrix[j, i] <- edge  # Ensure symmetry
          }
     }
     
     # Final binary conversion
     adjacency_matrix <- (adjacency_matrix + t(adjacency_matrix)) / 2
     adjacency_matrix[adjacency_matrix > 0] <- 1  
     
     # Convert to igraph object
     graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", diag = FALSE)
     
     # Return results
     return(list(
          adjacency_matrix = adjacency_matrix,
          probability_matrix = P,
          communities = node_communities,  
          community_sizes = community_sizes,
          graph = graph
     ))
}

get_vec_adjacency <- function(Y) 
{
     I <- nrow(Y)
     
     zro <- (1:I)[ Y %*% rep(1, I) == 0 ]
     if (length(zro) > 0) {
          Y <- Y[-zro, -zro]
          Ylabs <- Ylabs[-zro]
     }
     
     Y[is.na(Y)] <- 0
     Y[Y != 0 ]  <- 1
     
     Yvec <- as.matrix(rep(NA, I*(I-1)/2))
     for (i in 1:(I-1))
          for (ii in (i+1):I)
               Yvec[get_k(i, ii, I)] <- Y[i, ii]
     
     return(Yvec)
}

# Function to compute posterior mean of FDR and FNR
compute_posterior_fdr_fnr <- function(Xi_chain, Xi_true) 
{
     num_samples <- nrow(Xi_chain)  # Number of posterior samples
     fdr_values  <- numeric(num_samples)
     fnr_values  <- numeric(num_samples)
     
     for (i in seq_len(num_samples)) {
          tmp <- compute_fdr_fnr(Xi_chain[i, ], Xi_true)  # Compute FDR & FNR for each sample
          fdr_values[i] <- tmp$FDR
          fnr_values[i] <- tmp$FNR
     }
     
     # Compute posterior mean of FDR and FNR
     posterior_fdr <- mean(fdr_values)
     posterior_fnr <- mean(fnr_values)
     
     return(list(posterior_FDR = posterior_fdr, posterior_FNR = posterior_fnr))
}