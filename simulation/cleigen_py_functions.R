get_hyperpars <- function() 
{
     
     a_tau <- 3
     b_tau <- 2

     a_sig <- 3
     b_sig <- 2

     del2_U      <- 0.1  
     n_U         <- 0    
     n_tun_U     <- 100  
     
     del2_zeta   <- 0.1  
     n_zeta      <- 0    
     n_tun_zeta  <- 100  
     
     list(
          a_tau = a_tau, b_tau = b_tau,  
          a_sig = a_sig, b_sig = b_sig,  
          del2_U = del2_U, n_U = n_U, n_tun_U = n_tun_U,  
          del2_zeta = del2_zeta, n_zeta = n_zeta, n_tun_zeta = n_tun_zeta
     )
}

get_initial_values <- function(I, K, Q, hyps) 
{
     sigma <- log(K) / log(I)
     alpha <- sigma * K / I^sigma 
     
     Xi <- matrix(c(rep(0, I %% K), rep(seq_len(K) - 1, each = floor(I / K))), ncol = 1)
     tausq <- 1 / rgamma(1, shape = hyps$a_tau, rate = hyps$b_tau)
     sigsq <- 1 / rgamma(1, shape = hyps$a_sig, rate = hyps$b_sig)
     zeta  <- rnorm(1, mean = 0, sd = sqrt(tausq))
     U <- matrix(rnorm(K * Q, mean = 0, sd = sqrt(sigsq)), K, Q)
     
     list(
          alpha = alpha,
          sigma = sigma,
          Xi = Xi, 
          tausq = tausq, 
          sigsq = sigsq,
          zeta = zeta, 
          U = U)
}

get_chains_data <- function(I, K, Q, n_sams) 
{
     
     Eta_chain    <- vector("list", n_sams)
     Xi_chain     <- matrix(NA, n_sams, I)
     zeta_chain   <- numeric(n_sams)
     U_chain      <- vector("list", n_sams)
     tausq_chain  <- numeric(n_sams)
     sigsq_chain  <- numeric(n_sams)
     loglik_chain <- numeric(n_sams)
     
     list(
          Eta_chain    = Eta_chain, 
          Xi_chain     = Xi_chain, 
          zeta_chain   = zeta_chain, 
          U_chain      = U_chain,
          tausq_chain  = tausq_chain, 
          sigsq_chain  = sigsq_chain,
          loglik_chain = loglik_chain
     )
}

MCMC <- function(Y, K, Q, n_sams, n_burn, n_skip) 
{
     I     <- get_I(Y)
     THETA <- get_chains_data(I, K, Q, n_sams)
     
     hyps    <- get_hyperpars()
     a_tau   <- hyps$a_tau
     b_tau   <- hyps$b_tau
     a_sig   <- hyps$a_sig
     b_sig   <- hyps$b_sig
     
     del2_zeta   <- hyps$del2_zeta
     n_zeta      <- hyps$n_zeta
     n_tun_zeta  <- hyps$n_tun_zeta
     del2_U      <- hyps$del2_U
     n_U         <- hyps$n_U
     n_tun_U     <- hyps$n_tun_U
     
     init   <- get_initial_values(I, K, Q, hyps)
     Xi     <- init$Xi
     zeta   <- init$zeta
     U      <- init$U
     tausq  <- init$tausq
     sigsq  <- init$sigsq
     alpha  <- init$alpha
     sigma  <- init$sigma    
     
     B      <- n_burn + n_skip * n_sams
     n_disp <- floor(0.1 * B)
     mean_K <- 0
     
     for (b in 1:B) {
          tmp        <- sample_zeta(b, n_tun_zeta, del2_zeta, n_zeta, n_burn, I, K, tausq, zeta, U, Xi, Y)
          zeta       <- tmp$zeta
          del2_zeta  <- tmp$del2_zeta
          n_zeta     <- tmp$n_zeta
          n_tun_zeta <- tmp$n_tun_zeta
          
          tmp        <- sample_U(b, n_tun_U, del2_U, n_U, n_burn, I, K, Q, sigsq, zeta, U, Xi, Y)
          U          <- tmp$U
          del2_U     <- tmp$del2_U
          n_U        <- tmp$n_U
          n_tun_U    <- tmp$n_tun_U
          
          # Update Xi and track cluster evolution
          Xi_out     <- sample_Xi(I, Q, alpha, sigma, zeta, sigsq, U, Xi, Y)
          Xi         <- Xi_out$Xi
          U          <- Xi_out$U
          K          <- Xi_out$K
          
          mean_K <- ((b - 1) * mean_K + K) / b
          
          tausq <- sample_tausq(a_tau, b_tau, zeta)
          sigsq <- sample_sigsq(K, Q, a_sig, b_sig, U)
          
          if ((b > n_burn) & (b %% n_skip == 0)) {
               Eta <- get_Eta(K, zeta, U)
               i <- (b - n_burn) / n_skip
               repeat { 
                    ll <- loglik(I, K, c(Eta), Xi, Y)
                    if(!is.na(ll)) {
                         break
                    }
               }
               THETA$loglik_chain[i]  <- ll
               THETA$Eta_chain[[i]]   <- c(Eta)
               THETA$Xi_chain[i, ]    <- c(Xi)
               THETA$zeta_chain[i]    <- c(zeta)
               THETA$U_chain[[i]]     <- U
               THETA$tausq_chain[i]   <- c(tausq)
               THETA$sigsq_chain[i]   <- c(sigsq)
          }
          if (b %% n_disp == 0) {
               cat(
                    "Progress: ", formatC(100 * b / B, digits = 1, format = "f"), "% complete\n",
                    "zeta  mixing rate = ", formatC(100 * n_zeta / b, digits = 2, format = "f"), "%, del2_zeta  = ", formatC(del2_zeta, digits = 5, format = "f"), "\n",
                    "U     mixing rate = ", formatC(100 * n_U / (b * mean_K), digits = 2, format = "f"), "%, del2_U     = ", formatC(del2_U, digits = 5, format = "f"), "\n",
                    "Current number of clusters K = ", K, "\n",
                    "---------------------------------------------------\n",
                    sep = ""
               )
          }
     }
     
     return(THETA) 
}

YPPP<- function(Yna, na_indices, K, Q, n_sams, n_burn, n_skip) 
{
     I     <- get_I(Y)
     THETA <- get_chains_data(I, K, Q, n_sams)
     
     hyps    <- get_hyperpars()
     a_tau   <- hyps$a_tau
     b_tau   <- hyps$b_tau
     a_sig   <- hyps$a_sig
     b_sig   <- hyps$b_sig
     
     del2_zeta   <- hyps$del2_zeta
     n_zeta      <- hyps$n_zeta
     n_tun_zeta  <- hyps$n_tun_zeta
     del2_U      <- hyps$del2_U
     n_U         <- hyps$n_U
     n_tun_U     <- hyps$n_tun_U
     
     init   <- get_initial_values(I, K, Q, hyps)
     Xi     <- init$Xi
     zeta   <- init$zeta
     U      <- init$U
     tausq  <- init$tausq
     sigsq  <- init$sigsq
     alpha  <- init$alpha
     sigma  <- init$sigma    
     
     B      <- n_burn + n_skip * n_sams
     n_disp <- floor(0.1 * B)
     
     # Posterior predictive probabilities
     y_ppp <- rep(0, sum(na_indices))
     
     for (b in 1:B) {
          # Sample parameters
          Yna_out <- sample_Y(I, c(get_Eta(K, zeta, U)), Xi, na_indices, Yna)
          Yna     <- Yna_out
          
          tmp        <- sample_zeta(b, n_tun_zeta, del2_zeta, n_zeta, n_burn, I, K, tausq, zeta, U, Xi, Y)
          zeta       <- tmp$zeta
          del2_zeta  <- tmp$del2_zeta
          n_zeta     <- tmp$n_zeta
          n_tun_zeta <- tmp$n_tun_zeta
          
          tmp        <- sample_U(b, n_tun_U, del2_U, n_U, n_burn, I, K, Q, sigsq, zeta, U, Xi, Y)
          U          <- tmp$U
          del2_U     <- tmp$del2_U
          n_U        <- tmp$n_U
          n_tun_U    <- tmp$n_tun_U

          # Update Xi and track cluster evolution
          Xi_out     <- sample_Xi(I, Q, alpha, sigma, zeta, sigsq, U, Xi, Y)
          Xi         <- Xi_out$Xi
          U          <- Xi_out$U
          K          <- Xi_out$K
          
          tausq <- sample_tausq(a_tau, b_tau, zeta)
          sigsq <- sample_sigsq(K, Q, a_sig, b_sig, U)
          
          # Posterior predictive probabilities
          if ((b > n_burn) & (b %% n_skip == 0)) {
               y_ppp <- y_ppp + Yna[na_indices] / n_sams
          }
     }
     
     return(y_ppp) 
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