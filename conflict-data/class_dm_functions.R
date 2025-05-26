get_hyperpars <- function() {
     # Hyperparameters
     mu_mu   <- 0
     sig2_mu <- 3
     a_sig   <- 3
     b_sig   <- 2
     
     # Metropolis-Hastings parameters
     del2_Lambda  <- 0.1
     n_Lambda     <- 0
     n_tun_Lambda <- 100
     
     # Return as a list
     list(
          mu_mu         = mu_mu,
          sig2_mu       = sig2_mu,
          a_sig         = a_sig,
          b_sig         = b_sig,
          del2_Lambda   = del2_Lambda,
          n_Lambda      = n_Lambda,
          n_tun_Lambda  = n_tun_Lambda
     )
}

get_initial_values <- function(I, K, hyps) {
     # Initialize
     alpha <- solve_alpha(K, I)
     omega <- rep(alpha / K, K)
     Xi <- matrix(c(rep(0, I %% K), rep(seq_len(K) - 1, each = floor(I / K))), ncol = 1)
     sigsq <- 1 / rgamma(1, shape = hyps$a_sig, rate = hyps$b_sig)
     mu <- rnorm(1, mean = hyps$mu_mu, sd = sqrt(hyps$sig2_mu))
     Lambda <- rnorm(K * (K + 1) / 2, mean = mu, sd = sqrt(sigsq))
     
     # Return as a list
     list(
          alpha  = alpha,
          omega  = omega,
          Xi     = Xi,
          sigsq  = sigsq,
          mu     = mu,
          Lambda = Lambda
     )
}

get_chains_data <- function(I, K, n_sams, n_burn, n_skip) {
     # Initialize chains as matrices filled with NA_real_
     Lambda_chain <- matrix(NA_real_, n_sams, K * (K + 1) / 2)
     Xi_chain     <- matrix(NA_real_, n_sams, I)
     loglik_chain <- matrix(NA_real_, n_burn + n_skip * n_sams, 1)
     # mu_chain     <- matrix(NA_real_, n_sams, 1)
     # sigsq_chain  <- matrix(NA_real_, n_sams, 1)
     # omega_chain  <- matrix(NA_real_, n_sams, K)

     # Return list of chains
     list(
          Lambda_chain = Lambda_chain,
          Xi_chain     = Xi_chain,
          loglik_chain = loglik_chain
          # mu_chain     = mu_chain,
          # sigsq_chain  = sigsq_chain,
          # omega_chain  = omega_chain,
     )
}

solve_alpha <- function(EK_target, I) {
     # Define the function to find the root
     f <- function(alpha) alpha * log((alpha + I) / alpha) - EK_target
     
     # Solve for alpha using uniroot
     alpha_solution <- uniroot(f, lower = 1e-5, upper = I, tol = 1e-6)$root
     
     return(alpha_solution)
}

MCMC <- function(Y, K, n_sams, n_burn, n_skip) {
     # Get data dimensions and initialize storage
     I     <- get_I(Y)
     THETA <- get_chains_data(I, K, n_sams, n_burn, n_skip)
     
     # Retrieve hyperparameters and MH parameters
     hyps <- get_hyperpars()
     
     mu_mu   <- hyps$mu_mu
     sig2_mu <- hyps$sig2_mu
     a_sig   <- hyps$a_sig
     b_sig   <- hyps$b_sig
     
     del2_Lambda  <- hyps$del2_Lambda
     n_Lambda     <- hyps$n_Lambda
     n_tun_Lambda <- hyps$n_tun_Lambda
     
     # Initialize values
     init   <- get_initial_values(I, K, hyps)
     
     Lambda <- init$Lambda
     mu     <- init$mu
     sigsq  <- init$sigsq
     Xi     <- init$Xi
     omega  <- init$omega
     alpha  <- init$alpha
     
     # Define total iterations and display frequency
     B      <- n_burn + n_skip * n_sams
     n_disp <- floor(0.1 * B)
     
     # MCMC Sampling
     for (b in seq_len(B)) {
          # Sample parameters
          Lambda_out   <- sample_Lambda(b, n_tun_Lambda, del2_Lambda, n_Lambda, n_burn, I, K, sigsq, mu, Lambda, Xi, Y)
          Lambda       <- Lambda_out$Lambda
          del2_Lambda  <- Lambda_out$del2_Lambda
          n_Lambda     <- Lambda_out$n_Lambda
          n_tun_Lambda <- Lambda_out$n_tun_Lambda
          
          Xi    <- sample_Xi(I, K, omega, Lambda, Xi, Y)
          omega <- sample_omega(K, alpha, Xi)
          
          mu    <- sample_mu(K, mu_mu, sig2_mu, sigsq, Lambda)
          sigsq <- sample_sigsq(K, a_sig, b_sig, mu, Lambda)
          
          THETA$loglik_chain[b] <- loglik(I, K, Lambda, Xi, Y)

          # Store sampled values
          if (b > n_burn && (b %% n_skip == 0)) {
               i <- (b - n_burn) / n_skip
               THETA$Lambda_chain[i, ] <- Lambda
               THETA$Xi_chain[i, ]     <- Xi
               # THETA$mu_chain[i]       <- mu
               # THETA$sigsq_chain[i]    <- sigsq
               # THETA$omega_chain[i, ]  <- omega
          }
          
          # Display progress
          if (b %% n_disp == 0) {
               cat(
                    "Progress: ", formatC(100 * b / B, digits = 1, format = "f"), "% complete\n",
                    "Lambda mixing rate = ", formatC(100 * n_Lambda / (b * 0.5 * K * (K + 1)), digits = 2, format = "f"), "%, del2_Lambda = ", formatC(del2_Lambda, digits = 5, format = "f"), "\n",
                    "---------------------------------------------------\n",
                    sep = ""
               )
          }
     }
     
     return(THETA)
}

incidence_matrix <- function(THETA) {
     # Extract dimensions
     I <- ncol(THETA$Xi_chain)
     B <- nrow(THETA$Xi_chain)
     
     # Compute incidence matrix vector
     A_vec <- incidence_matrix0(I, B, THETA$Xi_chain)
     
     # Initialize symmetric matrix
     A <- matrix(0, I, I)
     
     # Fill upper triangle and copy to lower triangle
     for (i in seq_len(I - 1)) {
          for (ii in seq(i + 1, I)) {
               k <- get_k(i, ii, I)
               A[i, ii] <- A_vec[k]
               A[ii, i] <- A_vec[k]
          }
     }
     
     # Set diagonal elements to 1
     diag(A) <- 1
     
     return(A)
}