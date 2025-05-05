#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;

uword get_k(const int& i, const int& ii, const double& I) 
{
     // Precompute reusable terms
     double term1 = I * (I - 1.0) / 2.0;  // Total number of pairs
     double term2 = (I - i) * (I - i - 1.0) / 2.0;  // Remaining pairs after row i
     return static_cast<uword>(term1 - term2 + ii - i - 1.0);
}

uword get_k_diag(const int& k, const int& kk, const double& K) 
{
     // Precompute reusable terms
     double term1 = K * k;  // Base index shift
     double term2 = 0.5 * k * (k + 1.0);  // Triangular number adjustment
     return static_cast<uword>(term1 - term2 + kk);
}

void tunning(double& del2, int& n_tun, const double& mix_rate, const int& b) 
{
     // Declare constants within the function
     const double target_rate = 0.35;
     const double tolerance = 0.05;
     const double max_del2 = 10.0;
     const double min_del2 = 1e-10;
     
     if (b % n_tun == 0) {
          double diff = mix_rate - target_rate;
          if (std::abs(diff) > tolerance) {
               del2 *= std::exp(0.1 * diff);
               del2 = std::clamp(del2, min_del2, max_del2);
               n_tun = 100;
          } else {
               n_tun += 100;
          }
     }
}

uword wsample(vec probs) 
{
     // Normalize probabilities
     probs /= accu(probs);
     
     // Generate cumulative sum
     vec probsum = cumsum(probs);
     
     // Draw a uniform random number
     double u = R::runif(0.0, 1.0);
     
     // Use binary search for efficient lookup (O(log n))
     return std::lower_bound(probsum.begin(), probsum.end(), u) - probsum.begin();
}

vec rdirichlet(const vec& alpha) 
{
     uword K = alpha.n_elem;
     vec out(K);
     
     for (uword k = 0; k < K; k++) {
          out[k] = R::rgamma(alpha[k], 1.0);
     }
     
     return out / accu(out);
}

// [[Rcpp::export]]
double loglik(const double& I, const double& K, const vec& Eta, const uvec& Xi, const vec& Y) 
{
     double out = 0.0, prob;
     uword eta_index, y_index;
     
     // Iterate over all node pairs (i, j) to compute likelihood
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               // Compute indices
               eta_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
               y_index = get_k(i, ii, I);
               
               // Compute probability using the probit function
               prob = R::pnorm(Eta[eta_index], 0.0, 1.0, 1, 0);
               
               // Accumulate log-likelihood
               out += R::dbinom(Y[y_index], 1, prob, 1);
          }
     }
     
     return out;
}

// [[Rcpp::export]]
vec get_Eta(const double& K, const double& zeta, const mat& U) 
{
     vec out(0.5 * K * (K + 1.0));
     uword eta_index;
     
     // Compute eta_{k,l} using optimized eigenmodel computation
     for (uword k = 0; k < K; k++) {
          for (uword l = k; l < K; l++) {
               eta_index = get_k_diag(k, l, K);
               out[eta_index] = zeta + sum(U.row(k) % U.row(l));
          }
     }
     
     return out;
}

// [[Rcpp::export]]
List sample_zeta(const double& b, int n_tun_zeta, double del2_zeta, int n_zeta, 
                 const int& n_burn, const double& I, const double& K, const double& tausq, 
                 double zeta, const mat& U, const uvec& Xi, const vec& Y) 
{
     // Propose new zeta using Gaussian random walk
     double zeta_p = R::rnorm(zeta, sqrt(del2_zeta));
     
     // Compute log-likelihoods for current and proposed values
     double loglik_current = loglik(I, K, get_Eta(K, zeta, U), Xi, Y);
     double loglik_proposed = loglik(I, K, get_Eta(K, zeta_p, U), Xi, Y);
     
     // Compute log acceptance probability
     double log_accept_prob = loglik_proposed - loglik_current - (pow(zeta_p, 2) - pow(zeta, 2)) / (2.0 * tausq);
     
     // Accept or reject the proposed update
     if (R::runif(0, 1) < std::exp(log_accept_prob)) {
          zeta = zeta_p;
          n_zeta++;
     }
     
     // Adaptive tuning of proposal variance during burn-in phase
     if (b < n_burn) {
          double mix_rate = static_cast<double>(n_zeta) / b;
          tunning(del2_zeta, n_tun_zeta, mix_rate, b);
     }
     
     return List::create(Named("zeta") = zeta,
                         Named("del2_zeta") = del2_zeta,
                         Named("n_zeta") = n_zeta,
                         Named("n_tun_zeta") = n_tun_zeta);
}

double lfcd_U(const rowvec& x, const uword& k, const double& I, const double& K, 
              const double& sigsq, const double& zeta, mat U,
              const uvec& Xi, const vec& Y) 
{
     double out = -pow(norm(x, 2), 2) / (2.0 * sigsq);  // Prior term for u_k
     
     // Update latent position
     U.row(k) = x;
     vec Eta = get_Eta(K, zeta, U);
     uword eta_index, y_index;
     
     // Compute likelihood contribution
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               if ((Xi[i] == k) || (Xi[ii] == k)) {
                    eta_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
                    y_index = get_k(i, ii, I);
                    
                    // Directly use R::pnorm(...) in log-likelihood computation
                    out += R::dbinom(Y[y_index], 1, R::pnorm(Eta[eta_index], 0.0, 1.0, 1, 0), 1);
               }
          }
     }
     
     return out;
}

// [[Rcpp::export]]
List sample_U(const double& b, int n_tun_U, double del2_U, int n_U, const int& n_burn, 
              const double& I, const double& K, const double& Q, const double& sigsq, 
              const double& zeta, mat U, const uvec& Xi, const vec& Y) 
{
     rowvec u_c(Q), u_p(Q);
     double lfcd_current, lfcd_proposed;
     
     // Metropolis-Hastings step for each cluster k
     for (uword k = 0; k < K; k++) {
          u_c = U.row(k);
          u_p = u_c + sqrt(del2_U) * randn<rowvec>(Q);  // Propose new latent position
          
          // Compute log full conditional densities
          lfcd_current = lfcd_U(u_c, k, I, K, sigsq, zeta, U, Xi, Y);
          lfcd_proposed = lfcd_U(u_p, k, I, K, sigsq, zeta, U, Xi, Y);
          
          // Accept or reject the proposed update
          if (R::runif(0, 1) < std::exp(lfcd_proposed - lfcd_current)) {
               U.row(k) = u_p;
               n_U++;
          }
     }
     
     // Adaptive tuning of proposal variance during burn-in phase
     if (b < n_burn) {
          double mix_rate = static_cast<double>(n_U) / (b * K);
          tunning(del2_U, n_tun_U, mix_rate, b);
     }
     
     return List::create(Named("U") = U,
                         Named("del2_U") = del2_U,
                         Named("n_U") = n_U,
                         Named("n_tun_U") = n_tun_U);
}

// [[Rcpp::export]]
vec sample_omega(const double& K, const double& alpha, const uvec& Xi) {
     vec Alpha(K, fill::value(alpha / K));
     uvec id;
     
     // Count occurrences of each cluster and update Alpha
     for (uword k = 0; k < K; k++) {
          id = find(Xi == k);
          Alpha[k] += id.n_elem;
     }
     
     // Sample from Dirichlet distribution
     return rdirichlet(Alpha);
}

// [[Rcpp::export]]
uvec sample_Xi(const double& I, const double& K, const vec& omega, const vec& Eta, uvec Xi, const vec& Y) {
     uword y_index, eta_index;
     
     for (uword i = 0; i < I; i++) {
          vec logprobs(K, fill::zeros);  // Initialize log-probabilities
          
          for (uword k = 0; k < K; k++) {
               logprobs[k] = log(omega[k]);  // Log-prior contribution
               
               // Compute likelihood contributions
               if (i < I - 1) {
                    for (uword ii = i + 1; ii < I; ii++) {
                         y_index = get_k(i, ii, I);
                         eta_index = get_k_diag(std::min(k, Xi[ii]), std::max(k, Xi[ii]), K);
                         logprobs[k] += R::dbinom(Y[y_index], 1, R::pnorm(Eta[eta_index], 0.0, 1.0, 1, 0), 1);
                    }
               }
               
               if (i > 0) {
                    for (uword ii = 0; ii < i; ii++) {
                         y_index = get_k(ii, i, I);
                         eta_index = get_k_diag(std::min(k, Xi[ii]), std::max(k, Xi[ii]), K);
                         logprobs[k] += R::dbinom(Y[y_index], 1, R::pnorm(Eta[eta_index], 0.0, 1.0, 1, 0), 1);
                    }
               }
          }
          
          // Normalize and sample a new cluster assignment
          Xi[i] = wsample(exp(logprobs - logprobs.max()));
     }
     
     return Xi;
}

// [[Rcpp::export]]
double sample_sigsq(const double& K, const double& Q, const double& a_sig, 
                    const double& b_sig, const mat& U) 
{
     // Compute the shape and scale parameters for the Inverse-Gamma distribution
     double shape = a_sig + 0.5 * K * Q;
     double scale = 1.0 / (b_sig + 0.5 * accu(square(U)));  // Using square(U) for clarity
     
     // Sample from Inverse-Gamma distribution
     return 1.0 / R::rgamma(shape, scale);
}

// [[Rcpp::export]]
double sample_tausq(const double& a_tau, const double& b_tau, const double& zeta) 
{
     // Compute the shape and scale parameters for the Inverse-Gamma distribution
     double shape = a_tau + 0.5;
     double scale = 1.0 / (b_tau + 0.5 * pow(zeta, 2));
     
     // Sample from Inverse-Gamma distribution
     return 1.0 / R::rgamma(shape, scale);
}

// [[Rcpp::export]]
rowvec incidence_matrix0(const double& I, const double& B, const umat& Xi_chain) 
{
     rowvec out(0.5 * I * (I - 1.0), fill::zeros);
     const double scale_factor = 1.0 / B;
     uword k;
     
     for (uword b = 0; b < B; b++) {
          for (uword i = 0; i < I - 1; i++) {
               for (uword ii = i + 1; ii < I; ii++) {
                    if (Xi_chain(b, i) == Xi_chain(b, ii)) {
                         k = get_k(i, ii, I);
                         out(k) += scale_factor;
                    }
               }
          }
     }
     
     return out;
}

// [[Rcpp::export]]
rowvec interaction_probs0(const double& I, const double& K, const double& B, const mat& Lambda_chain, const umat& Xi_chain) {
     const uword N = 0.5 * I * (I - 1.0);
     rowvec out(N, fill::zeros);
     const double scale_factor = 1.0 / B;
     
     urowvec Xi(I);
     rowvec Lambda(0.5 * K * (K + 1.0));
     uword k, lambda_index;
     
     for (uword b = 0; b < B; b++) {
          Lambda = Lambda_chain.row(b);
          Xi = Xi_chain.row(b);
          
          for (uword i = 0; i < I - 1; i++) {
               for (uword ii = i + 1; ii < I; ii++) {
                    k = get_k(i, ii, I);
                    lambda_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
                    out(k) += R::pnorm(Lambda[lambda_index], 0.0, 1.0, 1, 0) * scale_factor;
               }
          }
     }
     return out;
}

// [[Rcpp::export]]
rowvec interaction_probs_iteration(const double& I, const double& K, const rowvec& Lambda, const urowvec& Xi) {
     const uword N = 0.5 * I * (I - 1.0);
     rowvec out(N, fill::zeros);
     
     uword k, lambda_index;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               k = get_k(i, ii, I);
               lambda_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
               out(k) +=R::pnorm(Lambda[lambda_index], 0.0, 1.0, 1, 0);
          }
     }
     
     return out;
}

// [[Rcpp::export]]
List WAIC(const double& I, const double& K, const double& B, 
          const vec& Y, const mat& Eta_chain, const umat& Xi_chain) 
{
     const double scale_factor = 1.0 / B;  // Precompute scaling factor
     
     uword m, eta_index;
     double tmp, a_ib, a_ib_ssq, a_ib_sum;
     
     urowvec Xi(I);
     rowvec Eta(0.5 * K * (K + 1.0));
     double lppd = 0.0, slp = 0.0, pWAIC2 = 0.0;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               m = get_k(i, ii, I);  // Index for Y
               tmp = 0.0;
               a_ib_ssq = 0.0;
               a_ib_sum = 0.0;
               
               for (uword b = 0; b < B; b++) {
                    Eta = Eta_chain.row(b);
                    Xi = Xi_chain.row(b);
                    
                    eta_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
                    double prob = R::pnorm(Eta[eta_index], 0.0, 1.0, 1, 0);
                    a_ib = R::dbinom(Y[m], 1, prob, 1);
                    
                    // WAIC 1: Sum of exponentiated log-likelihoods
                    tmp += std::exp(a_ib) * scale_factor;
                    slp += a_ib * scale_factor;
                    
                    // WAIC 2: Compute variance-related term
                    a_ib_ssq += pow(a_ib, 2);
                    a_ib_sum += a_ib;
               }
               
               // Compute log pointwise predictive density
               lppd += log(tmp);
               
               // Compute pWAIC2 (alternative measure of effective parameters)
               pWAIC2 += (a_ib_ssq - B * pow(a_ib_sum * scale_factor, 2)) / (B - 1.0);
          }
     }
     
     // Compute WAIC
     double pWAIC1 = 2.0 * lppd - 2.0 * slp;
     double waic1 = -2.0 * lppd + 2.0 * pWAIC1;
     double waic2 = -2.0 * lppd + 2.0 * pWAIC2;
     
     return List::create(Named("lppd")   = lppd,
                         Named("pWAIC1") = pWAIC1,
                         Named("pWAIC2") = pWAIC2,
                         Named("waic1")  = waic1,
                         Named("waic2")  = waic2);
}

// [[Rcpp::export]]
mat simulate_data(const double& I, const double& K, const vec& Eta, const uvec& Xi) {
     mat Y(I, I, fill::zeros);  // Initialize adjacency matrix
     uword eta_index;
     double prob;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               eta_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
               prob = R::pnorm(Eta[eta_index], 0.0, 1.0, 1, 0);
               
               Y(i, ii) = R::rbinom(1, prob);
               Y(ii, i) = Y(i, ii);  // Ensure symmetry
          }
     }
     
     return Y;
}

// [[Rcpp::export]]
arma::mat co_clustering_matrix(const arma::uvec& partition)
{
     uword n = partition.n_elem;
     arma::mat C = arma::zeros<arma::mat>(n, n);  // Initialize matrix with zeros
     
     // Construct co-clustering matrix
     for (uword i = 0; i < n; i++) {
          for (uword j = i; j < n; j++) {
               if (partition[i] == partition[j]) {
                    C(i, j) = 1;
                    C(j, i) = 1;  // Ensure symmetry
               }
          }
     }
     
     return C;
}

// [[Rcpp::export]]
List compute_fdr_fnr(const arma::uvec& partition_Z, const arma::uvec& partition_Z0) 
{
     // Compute co-clustering matrices
     arma::mat C_Z  = co_clustering_matrix(partition_Z);
     arma::mat C_Z0 = co_clustering_matrix(partition_Z0);
     
     // Compute True Positives (TP), False Positives (FP), False Negatives (FN)
     double TP = 0.5 * accu((C_Z == 1) % (C_Z0 == 1));  // True Positives
     double FP = 0.5 * accu((C_Z == 1) % (C_Z0 == 0));  // False Positives
     double FN = 0.5 * accu((C_Z == 0) % (C_Z0 == 1));  // False Negatives
     
     // Compute FDR and FNR
     double FDR = (TP + FP > 0) ? (FP / (TP + FP)) : 0;
     double FNR = (TP + FN > 0) ? (FN / (TP + FN)) : 0;
     
     return List::create(Named("FDR") = FDR, Named("FNR") = FNR);
}