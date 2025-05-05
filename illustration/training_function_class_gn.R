run_training <- function(route, data, model, prior, n_sams, n_burn, n_skip) {
     
     # ------------------------------------------------------------
     # Load required libraries
     # ------------------------------------------------------------
     suppressMessages(suppressWarnings({
          library(Rcpp)
          library(igraph)
          library(coda)
          library(corrplot)
          library(mcclust.ext)  # devtools::install_github("muschellij2/mcclust.ext")
     }))
     
     # ------------------------------------------------------------
     # Load data and source C++/R functions
     # ------------------------------------------------------------
     load(file = paste0(route, data, ".RData"))
     
     sourceCpp(paste0(route, "class_gn_functions.cpp"))
     
     r_files <- c(
          "r_functions.R", 
          "class_gn_functions.R", 
          "test_statistics_functions.R", 
          "esbm.R"
     )
     lapply(r_files, source)
     
     # ------------------------------------------------------------
     # Prepare file prefix and format adjacency matrix
     # ------------------------------------------------------------
     file_prefix <- paste0(data, "_", model, "_", prior)
     
     # ------------------------------------------------------------
     # Compute basic network characteristics
     # ------------------------------------------------------------
     node_degrees <- rowSums(Ycube)
     avg_degree   <- mean(node_degrees)
     K <- floor(I / avg_degree)  # Set number of communities
     
     # ------------------------------------------------------------
     # Run MCMC
     # ------------------------------------------------------------
     cat("Starting MCMC sampling...\n")
     set.seed(42)
     samples <- MCMC(Y, K, n_sams, n_burn, n_skip)
     cat("MCMC sampling completed.\n")
     
     # ------------------------------------------------------------
     # Plot log-likelihood chain
     # ------------------------------------------------------------
     loglik_clean <- samples$loglik_chain[
          !is.na(samples$loglik_chain) & samples$loglik_chain != -Inf
     ]
     
     pdf(file = paste0(route, file_prefix, "_loglik_chain.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     plot(
          x = seq_along(loglik_clean), y = loglik_clean, 
          type = "p", pch = ".", 
          ylim = quantile(loglik_clean, c(0.005, 0.995)), 
          xlab = "Iteration", ylab = "Log-likelihood"
     )
     dev.off()
     
     # ------------------------------------------------------------
     # Compute and plot interaction probabilities
     # ------------------------------------------------------------
     interaction_probs <- incidence_matrix(samples)
     interaction_probs[interaction_probs < 0] <- 0
     interaction_probs[interaction_probs > 1] <- 1
     
     Xi_hat_VI <- minVI(interaction_probs)
     cluster_assignments <- Xi_hat_VI$cl
     ordering <- order(cluster_assignments)
     interaction_probs <- interaction_probs[ordering, ordering]
     
     pdf(file = paste0(route, file_prefix, "_interaction_probs.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     corrplot(
          corr = interaction_probs, is.corr = FALSE, 
          addgrid.col = NA, method = "color", 
          tl.pos = "n", cl.pos = "n",
          col = colorRampPalette(c("white", rev(heat.colors(100))))(100)
     )
     dev.off()
     
     # ------------------------------------------------------------
     # Compute test statistics
     # ------------------------------------------------------------
     cat("Computing test statistics...\n")
     set.seed(42)
     test_stats <- class_gn_test(Ycube, I, samples)
     cat("Test statistics computation completed.\n")
     
     # ------------------------------------------------------------
     # Save results
     # ------------------------------------------------------------
     save(
          samples,
          interaction_probs,
          test_stats,
          file = paste0(route, file_prefix, ".RData")
     )
}