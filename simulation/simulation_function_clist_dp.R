run_simulation <- function(model, prior, density_type, cluster_type, size, 
                           n_sams, n_burn, n_skip) {
     
     suppressMessages({
          suppressWarnings({
               library(Rcpp)
               library(igraph)
               library(coda)
               library(corrplot)
               library(mcclust.ext)  # devtools::install_github("muschellij2/mcclust.ext")
          })
     })

     # Load compiled C++ functions
     sourceCpp("clist_dp_functions.cpp")
     
     # Load R functions
     r_files <- c(
          "r_functions.R", 
          "clist_dp_functions.R", 
          "test_statistics_functions.R", 
          "esbm.R"
     )
     lapply(r_files, source)
     
     # --- Node and Community Settings ---
     if (size == "small") {
          num_nodes       <- 50
          num_communities <- 3     
     } else if (size == "medium") {
          num_nodes       <- 100
          num_communities <- 5     
     } else if (size == "large") {
          num_nodes       <- 200
          num_communities <- 10     
     } else {
          stop("Invalid size parameter. Use 'small', 'medium', or 'large'.")
     }
     
     # --- Density Settings ---
     if (density_type == "hom") {
          p_within_values  <- seq(0.50, 0.60, by = 0.02)
          p_between_values <- seq(0.05, 0.10, by = 0.01)     
     } else if (density_type == "het") {
          p_within_values  <- seq(0.50, 0.60, by = 0.02)
          p_between_values <- seq(0.10, 0.15, by = 0.01)     
     } else {
          stop("Invalid density_type. Use 'hom' or 'het'.")
     }
     
     # --- File Name Construction ---
     file_prefix <- paste0("simulation_", model, "_", prior, "_", density_type, "_", cluster_type, "_", size)
     
     # --- Generate Data ---
     set.seed(42)
     sbm_simulation <- generate_SBM_network(
          num_nodes, num_communities,
          p_within_values, p_between_values, 
          cluster_type
     )
     
     # --- Convert adjacency matrix to vector format ---
     Y <- get_vec_adjacency(Y = sbm_simulation$adjacency_matrix)
     I <- get_I(Y)
     
     # --- Compute node degrees and average degree ---
     node_degrees <- rowSums(sbm_simulation$adjacency_matrix)
     avg_degree   <- mean(node_degrees)
     
     # --- Model Fitting Settings ---
     Q <- 4
     K <- floor(I / avg_degree)
     
     # --- Run MCMC Sampling ---
     cat("Starting MCMC sampling...\n")
     set.seed(42)
     samples <- MCMC(Y, K, Q, n_sams, n_burn, n_skip)
     cat("MCMC sampling completed.\n")
     
     # --- Save Log-likelihood Chain Plot ---
     loglik_clean <- samples$loglik_chain[!is.na(samples$loglik_chain) & samples$loglik_chain != -Inf]
     pdf(file = paste0(file_prefix, "_loglik_chain.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     plot(
          x = seq_len(length(loglik_clean)), y = loglik_clean, 
          type = "p", pch = ".",
          ylim = quantile(loglik_clean, c(0.005, 0.995)), 
          xlab = "Iteration", ylab = "Log-likelihood", main = ""
     )
     dev.off()
     
     # --- Compute Interaction Probabilities ---
     interaction_probs <- incidence_matrix(samples) 
     
     # --- Save Interaction Probabilities Matrix Plot ---
     pdf(file = paste0(file_prefix, "_interaction_probs.pdf"))
     par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
     colorscale <- c("white", rev(heat.colors(100)))
     corrplot(
          corr = interaction_probs, is.corr = FALSE, 
          addgrid.col = NA, method = "color", 
          tl.pos = "n",  # remove text labels
          col = colorRampPalette(colorscale)(100),
          cl.pos = "n"
     )
     dev.off() 

     # --- Compute Test Statistics ---
     cat("Computing test statistics...\n")
     set.seed(42)
     test_stats <- clist_dp_test(Ycube = sbm_simulation$adjacency_matrix, I, samples) 
     cat("Test statistics computation completed.\n")
     
     # --- Save Results ---
     save(
          samples,
          interaction_probs,
          test_stats,
          file = paste0(file_prefix, ".RData")
     )
}