incidence_matrix <- function(THETA) {
     # Extract dimensions
     I <- ncol(THETA$Xi_chain)
     B <- nrow(THETA$Xi_chain)
     
     # Compute incidence vector
     A_vec <- incidence_matrix0(I, B, THETA$Xi_chain)
     
     # Initialize symmetric matrix
     A <- matrix(0, I, I)
     
     # Fill upper triangle and mirror to lower triangle
     for (i in seq_len(I - 1)) {
          for (j in seq(i + 1, I)) {
               k <- get_k(i, j, I)
               A[i, j] <- A_vec[k]
               A[j, i] <- A_vec[k]
          }
     }
     
     # Set diagonal to 1
     diag(A) <- 1
     
     return(A)
}

interaction_matrix <- function(THETA, K) {
     # Extract dimensions
     I <- ncol(THETA$Xi_chain)
     B <- nrow(THETA$Xi_chain)
     
     # Compute interaction probabilities vector
     A_vec <- interaction_probs0(I, K, B, THETA$Lambda_chain, THETA$Xi_chain)
     
     # Initialize symmetric matrix
     A <- matrix(0, I, I)
     
     # Fill upper triangle and mirror to lower triangle
     for (i in seq_len(I - 1)) {
          for (j in seq(i + 1, I)) {
               k <- get_k(i, j, I)
               A[i, j] <- A_vec[k]
               A[j, i] <- A_vec[k]
          }
     }
     
     return(A)
}

interaction_matrix_np <- function(THETA) {
     # Extract dimensions
     I <- ncol(THETA$Xi_chain)
     B <- nrow(THETA$Xi_chain)
     
     # Compute interaction probabilities vector
     A_vec <- interaction_probs0(I, B, THETA$Lambda_chain, THETA$Xi_chain)
     
     # Initialize symmetric matrix
     A <- matrix(0, I, I)
     
     # Fill upper triangle and mirror to lower triangle
     for (i in seq_len(I - 1)) {
          for (j in seq(i + 1, I)) {
               k <- get_k(i, j, I)
               A[i, j] <- A_vec[k]
               A[j, i] <- A_vec[k]
          }
     }
     
     return(A)
}

ordering <- function(AM, IM, PM, labels) {
     
     # Remove isolated nodes
     isolated_nodes <- colSums(AM) == 0
     
     AM <- AM[!isolated_nodes, !isolated_nodes]
     IM <- IM[!isolated_nodes, !isolated_nodes]
     PM <- PM[!isolated_nodes, !isolated_nodes]
     labels <- labels[!isolated_nodes]
     
     I <- length(labels)
     
     # Relabel clusters with consecutive integers
     labels_new <- rep(NA, I)
     labels_unique <- sort(unique(labels))
     
     count <- 1
     for (i in labels_unique) {
          labels_new[labels == i] <- count
          count <- count + 1
     }
     
     labels <- labels_new
     
     # Relabel clusters by decreasing size
     labels_new <- rep(NA, I)
     
     cluster_sizes <- table(labels)
     labels_unique <- as.numeric(names(sort(cluster_sizes, decreasing = TRUE)))
     
     count <- 1
     for (i in labels_unique) {
          labels_new[labels == i] <- count
          count <- count + 1
     }
     
     labels <- labels_new
     
     # Reorder matrices based on new labels
     ordered_nodes <- order(labels)
     
     return(list(
          AM = AM[ordered_nodes, ordered_nodes],
          IM = IM[ordered_nodes, ordered_nodes],
          PM = PM[ordered_nodes, ordered_nodes],
          labels = labels
     ))
}

plot_adjacency <- function(Y, labels) {
     
     # Plot adjacency matrix
     corrplot(Y,
              method = "color",
              col = colorRampPalette(c("white", "red"))(10),
              tl.pos = "n",    # Hide text labels
              cl.pos = "n",    # Hide color legend
              is.corr = FALSE) # Not a correlation matrix
     
     # Add cluster boundaries
     cluster_sizes <- table(labels)
     n_nodes <- length(labels)
     
     cuts_v <- c(0, cumsum(cluster_sizes)) + 0.5
     cuts_h <- c(0, cumsum(rev(cluster_sizes))) + 0.5
     
     for (v in cuts_v) {
          segments(x0 = v, y0 = 0.5, x1 = v, y1 = n_nodes + 0.5,
                   col = "gray50", lwd = 1)
     }
     
     for (h in cuts_h) {
          segments(x0 = 0.5, y0 = h, x1 = n_nodes + 0.5, y1 = h,
                   col = "gray50", lwd = 1)
     }
}

plot_incidence_probs <- function(P, labels) {
     
     # Plot interaction probability matrix
     corrplot(P,
              method = "color",
              col = c("white", rev(heat.colors(100))),
              tl.pos = "n",     # Hide text labels
              cl.pos = "r",     # Show color legend
              cl.length = 11,   # Number of labels in color legend
              is.corr = FALSE)  # Not a correlation matrix
     
     # Add cluster boundaries
     cluster_sizes <- table(labels)
     n_nodes <- length(labels)
     
     cuts_v <- c(0, cumsum(cluster_sizes)) + 0.5
     cuts_h <- c(0, cumsum(rev(cluster_sizes))) + 0.5
     
     for (v in cuts_v) {
          segments(x0 = v, y0 = 0.5, x1 = v, y1 = n_nodes + 0.5,
                   col = "gray50", lwd = 1)
     }
     
     for (h in cuts_h) {
          segments(x0 = 0.5, y0 = h, x1 = n_nodes + 0.5, y1 = h,
                   col = "gray50", lwd = 1)
     }
}

get_posterior_nclusters <- function(Xi_chain) {
     # Define number of iterations
     B <- nrow(Xi_chain)
     
     # Initialize vector for number of clusters
     n_clusters <- integer(B)
     
     # Loop over iterations
     for (b in seq_len(B)) {
          clustering <- Xi_chain[b, ]
          n_clusters[b] <- length(unique(clustering))
     }
     
     return(n_clusters)
}

pred_metrics <- function(Ycube, I, K, samples, base_path, model, prior) {
     
     # Load libraries and compiled C++ functions
     Rcpp::sourceCpp(paste0(base_path, model, "_", prior, "_", "functions.cpp"))
     require(pROC)
     
     # Preallocate matrix for test statistics
     B <- nrow(samples$Xi_chain)
     n_pred_stats <- 3
     pred_stats <- matrix(NA, B, n_pred_stats)
     
     # Simulation loop
     for (b in seq_len(B)) {
          
          # Extract parameters
          Lambda <- samples$Lambda_chain[b, ]
          Xi     <- samples$Xi_chain[b, ]
          
          # Compute interaction probabilities
          A_vec <- interaction_probs_iteration(I, K, Lambda, Xi)
          A <- matrix(0, I, I)
          
          # Fill upper triangle and mirror to lower triangle
          for (i in seq_len(I - 1)) {
               for (j in seq(i + 1, I)) {
                    k <- get_k(i, j, I)
                    A[i, j] <- A_vec[k]
                    A[j, i] <- A_vec[k]
               }
          }
          
          # Clip probabilities to avoid log(0)
          eps <- 1e-6
          A <- pmin(pmax(A, eps), 1 - eps)
          
          # Extract upper triangle entries
          upper_idx <- which(upper.tri(Ycube))
          y_true <- Ycube[upper_idx]
          y_pred <- A[upper_idx]
          
          # Compute test statistics
          mse <- mean((y_true - y_pred)^2)
          logloss <- -mean(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred))
          roc_obj <- roc(response = y_true, predictor = y_pred, quiet = TRUE)
          auc_value <- auc(roc_obj)
          
          # Store results
          pred_stats[b, ] <- c(mse, logloss, auc_value)
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat(sprintf("Simulation progress: %.1f%% completed\n", 100 * b / B))
          }
     }
     
     # Compute posterior means
     mean_mse <- mean(pred_stats[, 1])
     mean_logloss <- mean(pred_stats[, 2])
     mean_auc <- mean(pred_stats[, 3])
     
     # Compute 95% credible intervals
     ci_mse <- quantile(pred_stats[, 1], probs = c(0.025, 0.975))
     ci_logloss <- quantile(pred_stats[, 2], probs = c(0.025, 0.975))
     ci_auc <- quantile(pred_stats[, 3], probs = c(0.025, 0.975))
     
     # Return results as a named list
     list(
          mean_mse = mean_mse,
          mean_logloss = mean_logloss,
          mean_auc = mean_auc,
          ci_mse = ci_mse,
          ci_logloss = ci_logloss,
          ci_auc = ci_auc
     )
}

pred_metrics_np <- function(Ycube, I, samples, base_path, model, prior) {
     
     # Load libraries and compiled C++ functions
     Rcpp::sourceCpp(paste0(base_path, model, "_", prior, "_", "functions.cpp"))
     require(pROC)
     
     # Preallocate matrix for test statistics
     B <- nrow(samples$Xi_chain)
     n_pred_stats <- 3
     pred_stats <- matrix(NA, B, n_pred_stats)
     
     # Simulation loop
     for (b in seq_len(B)) {
          
          # Extract parameters
          Lambda <- c(samples$Lambda_chain[[b]])
          Xi     <- c(samples$Xi_chain[b, ])
          K      <- length(table(Xi))
          
          # Compute interaction probabilities
          A_vec <- interaction_probs_iteration(I, K, Lambda, Xi)
          A <- matrix(0, I, I)
          
          # Fill upper triangle and mirror to lower triangle
          for (i in seq_len(I - 1)) {
               for (j in seq(i + 1, I)) {
                    k <- get_k(i, j, I)
                    A[i, j] <- A_vec[k]
                    A[j, i] <- A_vec[k]
               }
          }
          
          # Clip probabilities to avoid log(0)
          eps <- 1e-6
          A <- pmin(pmax(A, eps), 1 - eps)
          
          # Extract upper triangle entries
          upper_idx <- which(upper.tri(Ycube))
          y_true <- Ycube[upper_idx]
          y_pred <- A[upper_idx]
          
          # Compute test statistics
          mse <- mean((y_true - y_pred)^2)
          logloss <- -mean(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred))
          roc_obj <- roc(response = y_true, predictor = y_pred, quiet = TRUE)
          auc_value <- auc(roc_obj)
          
          # Store results
          pred_stats[b, ] <- c(mse, logloss, auc_value)
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat(sprintf("Simulation progress: %.1f%% completed\n", 100 * b / B))
          }
     }
     
     # Compute posterior means
     mean_mse <- mean(pred_stats[, 1])
     mean_logloss <- mean(pred_stats[, 2])
     mean_auc <- mean(pred_stats[, 3])
     
     # Compute 95% credible intervals
     ci_mse <- quantile(pred_stats[, 1], probs = c(0.025, 0.975))
     ci_logloss <- quantile(pred_stats[, 2], probs = c(0.025, 0.975))
     ci_auc <- quantile(pred_stats[, 3], probs = c(0.025, 0.975))
     
     # Return results as a named list
     list(
          mean_mse = mean_mse,
          mean_logloss = mean_logloss,
          mean_auc = mean_auc,
          ci_mse = ci_mse,
          ci_logloss = ci_logloss,
          ci_auc = ci_auc
     )
}

get_results_dm <- function(base_path, model, prior, dataset, n_thin){
     # Load data
     load(paste0(base_path, dataset, ".RData"))
     load(paste0(base_path, dataset, "_", model, "_", prior, ".RData"))
     
     # Load C++ and R helper functions
     Rcpp::sourceCpp(paste0(base_path, model, "_", prior, "_", "functions.cpp"))
     source(paste0(base_path, "r_functions.R"))
     source(paste0(base_path, "esbm.R"))
     
     # Prepare file prefix
     file_prefix <- paste0(dataset, "_", model, "_", prior)
     
     if (model %in% c("clist", "cleigen")) {
          names(samples)[names(samples) == "Eta_chain"] <- "Lambda_chain"
     }
     
     # Define basic dimensions
     B <- dim(samples$Xi_chain)[1]  # Number of MCMC iterations
     I <- dim(samples$Xi_chain)[2]  # Number of nodes
     
     # Estimate number of communities
     node_degrees <- rowSums(Ycube)
     avg_degree <- mean(node_degrees)
     K <- floor(I / avg_degree)
     
     # Thinning MCMC chains
     if (is.na(NA)) { n_thin <- (B/2)/1000 }
     indices_keep <- seq(from = B/2 + n_thin, to = B, by = n_thin)
     B <- length(indices_keep)
     
     samples$Xi_chain      <- samples$Xi_chain[indices_keep, ]
     samples$Lambda_chain  <- samples$Lambda_chain[indices_keep, ]
     test_stats$test_stats <- test_stats$test_stats[indices_keep, ]
     
     # Compute estimated interaction probabilities
     incidence_probs <- incidence_matrix(samples)
     
     # Truncate probabilities to [0, 1]
     incidence_probs[incidence_probs < 0] <- 0
     incidence_probs[incidence_probs > 1] <- 1
     
     # Find cluster assignments via VI minimization
     Xi_hat_VI <- minVI(incidence_probs)
     cluster_assignments <- Xi_hat_VI$cl
     
     # Estimate interaction matrix
     interaction_probs <- interaction_matrix(samples, K)
     
     # Compute predictive statistics
     cat("Computing predictive performance metrics:", "\n")
     pred_performance_stats <- pred_metrics(Ycube, I, K, samples, base_path, model, prior)
     
     # Reorder matrices based on clustering
     sorted_network_data <- ordering(Ycube, incidence_probs, interaction_probs, cluster_assignments)
     
     # Fix scale for visualization
     sorted_network_data$PM[1, 1] <- 1
     
     # Update number of nodes
     I <- dim(sorted_network_data$AM)[1]
     
     # Plot data
     pdf(file = paste0(base_path, file_prefix, "_adjacency_matrix.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     plot_adjacency(sorted_network_data$AM, sorted_network_data$labels)
     dev.off()
     
     # Plot co-clustering probabilities
     pdf(file = paste0(base_path, file_prefix, "_incidence_matrix.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     plot_incidence_probs(sorted_network_data$IM, sorted_network_data$labels)
     dev.off()
     
     # Plot interaction probabilities
     pdf(file = paste0(base_path, file_prefix, "_interaction_matrix.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     plot_incidence_probs(sorted_network_data$PM, sorted_network_data$labels)
     dev.off()
     
     # Posterior number of clusters
     post_nclusters <- get_posterior_nclusters(samples$Xi_chain)
     
     # Compute modularity
     g <- graph_from_adjacency_matrix(sorted_network_data$AM, mode = "undirected", diag = FALSE)
     modularity_score <- modularity(x = g, membership = sorted_network_data$labels)
     
     # Compute WAIC
     waic <- WAIC(I, K, B, get_vec_adjacency(sorted_network_data$AM), samples$Lambda_chain, samples$Xi_chain)
     
     # return
     list(
          B = B,
          I = I,
          test_stats = test_stats,
          incidence_probs = incidence_probs,
          cluster_assignments = cluster_assignments,
          interaction_probs = interaction_probs,
          pred_performance_stats = pred_performance_stats,
          sorted_network_data = sorted_network_data,
          post_nclusters = post_nclusters,
          modularity_score = modularity_score,
          waic = waic
     )
}

get_results_np <- function(base_path, model, prior, dataset, n_thin){
     # Load data
     load(paste0(base_path, dataset, ".RData"))
     load(paste0(base_path, dataset, "_", model, "_", prior, ".RData"))
     
     # Load C++ and R helper functions
     Rcpp::sourceCpp(paste0(base_path, model, "_", prior, "_", "functions.cpp"))
     source(paste0(base_path, "r_functions.R"))
     source(paste0(base_path, "esbm.R"))
     
     # Prepare file prefix
     file_prefix <- paste0(dataset, "_", model, "_", prior)
     
     if (model %in% c("clist", "cleigen")) {
          names(samples)[names(samples) == "Eta_chain"] <- "Lambda_chain"
     }
     
     # Define basic dimensions
     B <- dim(samples$Xi_chain)[1]  # Number of MCMC iterations
     I <- dim(samples$Xi_chain)[2]  # Number of nodes
     
     # Thinning MCMC chains
     if (is.na(NA)) { n_thin <- (B/2)/1000 }
     indices_keep <- seq(from = n_thin, to = B, by = n_thin)
     B <- length(indices_keep)
     
     samples$Xi_chain <- samples$Xi_chain[indices_keep, ]
     samples$Lambda_chain <- samples$Lambda_chain[indices_keep]
     test_stats$test_stats <- test_stats$test_stats[indices_keep, ]
     
     # Compute estimated interaction probabilities
     incidence_probs <- incidence_matrix(samples)
     
     # Truncate probabilities to [0, 1]
     incidence_probs[incidence_probs < 0] <- 0
     incidence_probs[incidence_probs > 1] <- 1
     
     # Find cluster assignments via VI minimization
     Xi_hat_VI <- minVI(incidence_probs)
     cluster_assignments <- Xi_hat_VI$cl
     
     # Estimate interaction matrix
     interaction_probs <- interaction_matrix_np(samples)
     
     # Compute predictive statistics
     cat("Computing predictive performance metrics:", "\n")
     pred_performance_stats <- pred_metrics_np(Ycube, I, samples, base_path, model, prior)
     
     # Reorder matrices based on clustering
     sorted_network_data <- ordering(Ycube, incidence_probs, interaction_probs, cluster_assignments)
     
     # Fix scale for visualization
     sorted_network_data$PM[1, 1] <- 1
     
     # Update number of nodes
     I <- dim(sorted_network_data$AM)[1]
     
     # Plot data
     pdf(file = paste0(base_path, file_prefix, "_adjacency_matrix.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     plot_adjacency(sorted_network_data$AM, sorted_network_data$labels)
     dev.off()
     
     # Plot co-clustering probabilities
     pdf(file = paste0(base_path, file_prefix, "_incidence_matrix.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     plot_incidence_probs(sorted_network_data$IM, sorted_network_data$labels)
     dev.off()
     
     # Plot interaction probabilities
     pdf(file = paste0(base_path, file_prefix, "_interaction_matrix.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     plot_incidence_probs(sorted_network_data$PM, sorted_network_data$labels)
     dev.off()
     
     # Posterior number of clusters
     post_nclusters <- get_posterior_nclusters(samples$Xi_chain)
     
     # Compute modularity
     g <- graph_from_adjacency_matrix(sorted_network_data$AM, mode = "undirected", diag = FALSE)
     modularity_score <- modularity(x = g, membership = sorted_network_data$labels)
     
     # Compute WAIC
     waic <- WAIC(I, B, get_vec_adjacency(sorted_network_data$AM), samples$Lambda_chain, samples$Xi_chain)
     
     # return
     list(
          B = B,
          I = I,
          test_stats = test_stats,
          incidence_probs = incidence_probs,
          cluster_assignments = cluster_assignments,
          interaction_probs = interaction_probs,
          pred_performance_stats = pred_performance_stats,
          sorted_network_data = sorted_network_data,
          post_nclusters = post_nclusters,
          modularity_score = modularity_score,
          waic = waic
     )
}