get_hyperpars <- function()
{
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

get_initial_values <- function(I, K, hyps) 
{
     gamma <- K / (K + log(I)) 

     base_size <- floor(I / K)
     extra_nodes <- I %% K
     Xi <- rep(seq_len(K) - 1, times = base_size)
     if (extra_nodes > 0)
          Xi <- c(Xi, seq_len(extra_nodes) - 1)
     
     sigsq  <- 1 / rgamma(1, shape = hyps$a_sig, rate = hyps$b_sig)
     mu     <- rnorm(1, mean = hyps$mu_mu, sd = sqrt(hyps$sig2_mu))
     Lambda <- rnorm(K * (K + 1) / 2, mean = mu, sd = sqrt(sigsq))
     
     list(gamma = gamma, Xi = Xi, sigsq = sigsq, mu = mu, Lambda = Lambda)
}

get_chains_data <- function(I, n_sams, n_burn, n_skip)
{
     Lambda_chain <- vector("list", n_sams)
     Xi_chain     <- matrix(NA_integer_, n_sams, I)
     loglik_chain <- numeric(n_burn + n_skip * n_sams)
     # mu_chain     <- numeric(n_sams)
     # sigsq_chain  <- numeric(n_sams)
     
     list(
          Lambda_chain = Lambda_chain,
          Xi_chain     = Xi_chain,
          loglik_chain = loglik_chain
          # mu_chain     = mu_chain,
          # sigsq_chain  = sigsq_chain,
     )
}

MCMC <- function(Y, K, n_sams, n_burn, n_skip) 
{
     # Get data dimensions and initialize storage
     I     <- get_I(Y)
     THETA <- get_chains_data(I, n_sams, n_burn, n_skip)
     
     # Retrieve hyperparameters and MH parameters
     hyps    <- get_hyperpars()
     mu_mu   <- hyps$mu_mu
     sig2_mu <- hyps$sig2_mu
     a_sig   <- hyps$a_sig
     b_sig   <- hyps$b_sig

     del2_Lambda  <- hyps$del2_Lambda
     n_Lambda     <- hyps$n_Lambda
     n_tun_Lambda <- hyps$n_tun_Lambda

     # Initialize values
     init    <- get_initial_values(I, K, hyps)
     Lambda  <- init$Lambda
     mu      <- init$mu
     sigsq   <- init$sigsq
     Xi      <- init$Xi
     gamma   <- init$gamma

     # Define total iterations and display frequency
     B      <- n_burn + n_skip * n_sams
     n_disp <- floor(0.1 * B)
     
     # Mean number of elements of Lambda
     mean_K <- 0
     
     # MCMC Sampling
     for (b in 1:B) {
          # Sample parameters
          Xi_out <- sample_Xi(I, gamma, mu, sigsq, Lambda, Xi, Y) 
          Xi     <- Xi_out$Xi
          K      <- Xi_out$K
          Lambda <- Xi_out$Lambda
          mean_K <- ((b - 1) * mean_K + K)/b
          
          Lambda_out   <- sample_Lambda(b, n_tun_Lambda, del2_Lambda, n_Lambda, n_burn, mean_K, I, K, sigsq, mu, Lambda, Xi, Y)
          Lambda       <- Lambda_out$Lambda
          del2_Lambda  <- Lambda_out$del2_Lambda
          n_Lambda     <- Lambda_out$n_Lambda
          n_tun_Lambda <- Lambda_out$n_tun_Lambda

          mu     <- sample_mu(K, mu_mu, sig2_mu, sigsq, Lambda)
          sigsq  <- sample_sigsq(K, a_sig, b_sig, mu, Lambda)
          
          THETA$loglik_chain[b] <- loglik(I, K, Lambda, Xi, Y)

          # Store sampled values
          if ((b > n_burn) & (b %% n_skip == 0)) {
               i <- (b - n_burn) / n_skip
               THETA$Xi_chain[i, ]     <- Xi
               THETA$Lambda_chain[[i]] <- Lambda
               # THETA$mu_chain[i]       <- mu
               # THETA$sigsq_chain[i]    <- sigsq
          }
          
          # Display progress
          if (b %% n_disp == 0) {
               cat(
                    "Progress: ", formatC(100 * b / B, digits = 1, format = "f"), "% complete\n",
                    "Lambda mixing rate = ", formatC(100 * n_Lambda / (b * 0.5 * mean_K * (mean_K + 1)), digits = 2, format = "f"), "%, del2_Lambda = ", formatC(del2_Lambda, digits = 5, format = "f"), "\n",
                    "Current number of clusters K = ", K, "\n",
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
     
     # Fill upper triangle using nested for loops
     for (i in 1:(I - 1)) {
          for (ii in (i + 1):I) {
               k <- get_k(i, ii, I)
               A[i, ii] <- A_vec[k]
               A[ii, i] <- A[i, ii]
          }
     }
     
     diag(A) <- 1
     
     return(A)
}