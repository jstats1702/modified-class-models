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

uword get_k(const int& i, const int& ii, const double& I) {
     // Precompute reusable terms
     double term1 = I * (I - 1.0) / 2.0;  // Total number of pairs
     double term2 = (I - i) * (I - i - 1.0) / 2.0;  // Remaining pairs after row i
     return static_cast<uword>(term1 - term2 + ii - i - 1.0);
}

uword get_k_diag(const int& k, const int& kk, const double& K) {
     // Precompute reusable terms
     double term1 = K * k;  // Base index shift
     double term2 = 0.5 * k * (k + 1.0);  // Triangular number adjustment
     return static_cast<uword>(term1 - term2 + kk);
}

void tunning(double& del2, int& n_tun, const double& mix_rate, const int& b) {
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

uword wsample(vec probs) {
     // Normalize probabilities
     probs /= accu(probs);
     
     // Generate cumulative sum
     vec probsum = cumsum(probs);
     
     // Draw a uniform random number
     double u = R::runif(0.0, 1.0);
     
     // Use binary search for efficient lookup (O(log n))
     return std::lower_bound(probsum.begin(), probsum.end(), u) - probsum.begin();
}

vec rdirichlet(const vec& alpha) {
     uword K = alpha.n_elem;
     vec out(K);
     
     for (uword k = 0; k < K; k++) {
          out[k] = R::rgamma(alpha[k], 1.0);
     }
     
     return out / accu(out);
}

// [[Rcpp::export]]
double loglik(const double& I, const double& K, const vec& Lambda, const uvec& Xi, const vec& Y) {
     // Xi is an I x 1 vector whose indices are zero-based
     double out = 0.0, prob;
     uword y_index, lambda_index;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               y_index = get_k(i, ii, I);  // Compute once per iteration
               lambda_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
               
               // Compute probability
               prob = R::pnorm(Lambda[lambda_index], 0.0, 1.0, 1, 0);
               
               // Accumulate log-likelihood
               out += R::dbinom(Y[y_index], 1, prob, 1);
          }
     }
     return out;
}

// [[Rcpp::export]]
double sample_sigsq(const double& K, const double& a_sig, const double& b_sig, const double& mu, const vec& Lambda) {
     const double shape = a_sig + K * (K + 1.0) / 4.0;
     const double scale = 1.0 / (b_sig + 0.5 * accu(pow(Lambda - mu, 2)));
     
     return 1.0 / R::rgamma(shape, scale);
}

// [[Rcpp::export]]
double sample_mu(const double& K, const double& mu_mu, const double& sig2_mu, const double& sigsq, const vec& Lambda) {
     const double v2 = 1.0 / (1.0 / sig2_mu + (0.5 * K * (K + 1.0)) / sigsq);
     const double m = v2 * (mu_mu / sig2_mu + accu(Lambda) / sigsq);
     
     return R::rnorm(m, sqrt(v2));
}

// SAMPLER FOR alpha IF NEEDED
/*
 // [[Rcpp::export]]
 List sample_alpha(const double& b, int n_tun_alpha, double del2_alpha, int n_alpha, const int& n_burn, 
 const double& K, const double& a_alpha, const double& b_alpha, double alpha, const vec& omega) {
 
 // Metropolis step
 const double sqrt_del2_alpha = sqrt(del2_alpha);
 double theta_c = log(alpha);
 double theta_p, alpha_p;
 
 // Propose new alpha, ensuring alpha_p >= 0.1
 do {
 theta_p = R::rnorm(theta_c, sqrt_del2_alpha);
 alpha_p = exp(theta_p);
 } while (alpha_p < 0.1);
 
 // Compute log acceptance probability
 double log_accept_prob = 
 (lgamma(alpha_p) - lgamma(alpha)) - 
 K * (lgamma(alpha_p / K) - lgamma(alpha / K)) + 
 (alpha_p - alpha) * accu(log(omega)) / K + 
 (a_alpha - 1.0) * (log(alpha_p) - log(alpha)) - 
 b_alpha * (alpha_p - alpha) + 
 (theta_p - theta_c);
 
 // Accept or reject
 if (R::runif(0.0, 1.0) < exp(log_accept_prob)) {
 alpha = alpha_p;
 n_alpha++;
 }
 
 // Adaptive tuning if within burn-in period
 if (b < n_burn) {
 double mix_rate = static_cast<double>(n_alpha) / b;
 tunning(del2_alpha, n_tun_alpha, mix_rate, b);
 }
 
 // Return updated values as a list
 return List::create(
 Named("alpha")       = alpha,
 Named("del2_alpha")  = del2_alpha,
 Named("n_alpha")     = n_alpha,
 Named("n_tun_alpha") = n_tun_alpha
 );
 }
 */

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
uvec sample_Xi(const double& I, const double& K, const vec& omega, const vec& Lambda, uvec Xi, const vec& Y) {
     // Xi contains uwords from 0 (class 1) to K-1 (class K)
     uword y_index, lambda_index;
     
     for (uword i = 0; i < I; i++) {
          vec logprobs(K, fill::zeros);
          
          for (uword k = 0; k < K; k++) {
               logprobs[k] = log(omega[k]);  // Initialize log-probabilities
               
               // Compute likelihood contributions
               if (i < I - 1) {
                    for (uword ii = i + 1; ii < I; ii++) {
                         y_index = get_k(i, ii, I);
                         lambda_index = get_k_diag(std::min(k, Xi[ii]), std::max(k, Xi[ii]), K);
                         logprobs[k] += R::dbinom(Y[y_index], 1, R::pnorm(Lambda[lambda_index], 0.0, 1.0, 1, 0), 1);
                    }
               }
               
               if (i > 0) {
                    for (uword ii = 0; ii < i; ii++) {
                         y_index = get_k(ii, i, I);
                         lambda_index = get_k_diag(std::min(k, Xi[ii]), std::max(k, Xi[ii]), K);
                         logprobs[k] += R::dbinom(Y[y_index], 1, R::pnorm(Lambda[lambda_index], 0.0, 1.0, 1, 0), 1);
                    }
               }
          }
          
          // Normalize and sample new Xi[i]
          Xi[i] = wsample(exp(logprobs - logprobs.max()));
     }
     
     return Xi;
}

// [[Rcpp::export]]
List sample_Lambda(const double& b, int n_tun_Lambda, double del2_Lambda, int n_Lambda, const int& n_burn, 
                   const double& I, const double& K, const double& sigsq, const double& mu, vec Lambda, const uvec& Xi, const vec& Y) {
     
     uword m;  
     double skl, nkl, lambda_c, lambda_p, log_accept_prob;
     uvec ind_k, ind_l;
     
     for (uword k = 0; k < K; k++) {
          for (uword l = k; l < K; l++) {
               // Compute sufficient statistics
               skl = 0.0;
               nkl = 0.0;
               
               ind_k = find(Xi == k);  // Get indices of nodes in cluster k
               ind_l = find(Xi == l);  // Get indices of nodes in cluster l
               
               for (uword i = 0; i < ind_k.n_elem; i++) {
                    for (uword ii = 0; ii < ind_l.n_elem; ii++) {
                         if (ind_k[i] < ind_l[ii]) {
                              skl += static_cast<double>(Y[get_k(ind_k[i], ind_l[ii], I)]);
                              nkl++;
                         }
                    }
               }
               
               // Metropolis step
               m = get_k_diag(k, l, K);
               lambda_c = Lambda[m];
               lambda_p = R::rnorm(lambda_c, sqrt(del2_Lambda));
               
               // Compute log-acceptance probability
               double log_pnorm_c = R::pnorm(lambda_c, 0.0, 1.0, 1, 1);
               double log_pnorm_p = R::pnorm(lambda_p, 0.0, 1.0, 1, 1);
               
               log_accept_prob = skl * (log_pnorm_p - log_pnorm_c) + 
                    (nkl - skl) * (log1p(-exp(log_pnorm_p)) - log1p(-exp(log_pnorm_c))) + 
                    (-0.5 / sigsq) * (pow(lambda_p - mu, 2) - pow(lambda_c - mu, 2));
               
               if (R::runif(0, 1) < exp(log_accept_prob)) {
                    Lambda[m] = lambda_p;
                    n_Lambda++;
               }
          }
     }
     
     // Adaptive tuning if within burn-in period
     if (b < n_burn) {
          double mix_rate = static_cast<double>(n_Lambda) / (b * 0.5 * K * (K + 1.0));
          tunning(del2_Lambda, n_tun_Lambda, mix_rate, b);
     }
     
     return List::create(
          Named("Lambda")       = Lambda,
          Named("del2_Lambda")  = del2_Lambda,
          Named("n_Lambda")     = n_Lambda,
          Named("n_tun_Lambda") = n_tun_Lambda
     );
}

// [[Rcpp::export]]
rowvec incidence_matrix0(const double& I, const double& B, const umat& Xi_chain) {
     const uword N = 0.5 * I * (I - 1.0);
     rowvec out(N, fill::zeros);
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
List WAIC(const double& I, const double& K, const double& B, const vec& Y, const mat& Lambda_chain, const umat& Xi_chain) {
     uword m, lambda_index;
     double tmp, a_ib, a_ib_ssq, a_ib_sum;
     const double scale_factor = 1.0 / B;
     
     urowvec Xi(I);
     rowvec Lambda(0.5 * K * (K + 1.0));
     double lppd = 0.0, slp = 0.0, pWAIC2 = 0.0;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               m = get_k(i, ii, I);
               tmp = 0.0;
               a_ib_ssq = 0.0;
               a_ib_sum = 0.0;
               
               for (uword b = 0; b < B; b++) {
                    Lambda = Lambda_chain.row(b);
                    Xi = Xi_chain.row(b);
                    
                    lambda_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
                    double prob = R::pnorm(Lambda[lambda_index], 0.0, 1.0, 1, 0);
                    a_ib = R::dbinom(Y[m], 1, prob, 1);
                    
                    // WAIC 1
                    tmp += exp(a_ib) * scale_factor;
                    slp += a_ib * scale_factor;
                    
                    // WAIC 2
                    a_ib_ssq += pow(a_ib, 2);
                    a_ib_sum += a_ib;
               }
               
               lppd += log(tmp);
               pWAIC2 += (a_ib_ssq - B * pow(a_ib_sum * scale_factor, 2)) / (B - 1.0);
          }
     }
     
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
mat simulate_data(const double& I, const double& K, const vec& Lambda, const uvec& Xi) {
     mat Y(I, I, fill::zeros);
     uword lambda_idx;
     double prob;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               lambda_idx = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
               prob = R::pnorm(Lambda[lambda_idx], 0.0, 1.0, 1, 0);
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