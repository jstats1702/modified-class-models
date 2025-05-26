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
     // vectorize upper triangular by rows
     // Precompute reusable terms
     double term1 = I * (I - 1.0) / 2.0;  // Total number of pairs
     double term2 = (I - i) * (I - i - 1.0) / 2.0;  // Remaining pairs after row i
     return static_cast<uword>(term1 - term2 + ii - i - 1.0);
}

uword get_k_diag(const int& k, const int& kk, const double& K) 
{
     // vectorize upper triangular by rows
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

// [[Rcpp::export]]
double loglik(const double& I, const double& K, const vec& Lambda, const uvec& Xi, const vec& Y)
{
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
double sample_sigsq(const double& K, const double& a_sig, const double& b_sig, const double& mu, const vec& Lambda) 
{
     const double shape = a_sig + K * (K + 1.0) / 4.0;
     const double scale = 1.0 / (b_sig + 0.5 * accu(pow(Lambda - mu, 2)));
     
     return 1.0 / R::rgamma(shape, scale);
}

// [[Rcpp::export]]
double sample_mu(const double& K, const double& mu_mu, const double& sig2_mu, const double& sigsq, const vec& Lambda) 
{
     const double v2 = 1.0 / (1.0 / sig2_mu + (0.5 * K * (K + 1.0)) / sigsq);
     const double m = v2 * (mu_mu / sig2_mu + accu(Lambda) / sigsq);
     
     return R::rnorm(m, sqrt(v2));
}


// [[Rcpp::export]]
List sample_Xi(const double& I, const double& gamma, const double& mu, 
                   const double& sigsq, vec Lambda, uvec Xi, const vec& Y) {
     
     uword y_index, lambda_index, k_old, kk_old, K = Xi.max() + 1;  // Initial number of clusters (zero-based)
     double log_prior_term;
     
     // Iterate over each node
     for (uword i = 0; i < I; i++) {
          
          // Compute cluster sizes
          arma::uvec cluster_sizes = arma::hist(Xi, K);  
          uword K_star = sum(cluster_sizes > 0);  // Number of non-empty clusters
          
          // Compute log probabilities for each cluster
          arma::vec logprobs(K + 1, fill::zeros);
          
          for (uword k = 0; k < K; k++) {
               if (cluster_sizes[k] > 0) {
                    logprobs[k] = log((cluster_sizes[k] + 1) * (I - K_star + gamma));
               } else {
                    logprobs[k] = -std::numeric_limits<double>::infinity(); // Avoid assigning to empty clusters
               }
               
               for (uword j = 0; j < I; j++) {
                    if (i == j) 
                         continue;
                    
                    y_index = (i < j) ? get_k(i, j, I) : get_k(j, i, I);
                    lambda_index = get_k_diag(std::min(k, Xi[j]), std::max(k, Xi[j]), K);
                    
                    logprobs[k] += R::dbinom(Y[y_index], 1, R::pnorm(Lambda[lambda_index], 0.0, 1.0, 1, 0), 1);
               }
          }
          
          // Compute probability for a new cluster under the Gnedin prior
          if (K_star > 1) {
               logprobs[K] = log(K_star * (K_star - gamma));
          } else {
               logprobs[K] = -std::numeric_limits<double>::infinity();  // Prevent forming a new cluster if only one exists
          }
          
          // Create a new Lambda matrix for the updated cluster assignments
          uword new_K = K + 1;
          arma::vec new_Lambda(new_K * (new_K + 1) / 2, fill::zeros);
          
          for (uword k = 0; k < new_K; k++) {
               for (uword kk = k; kk < new_K; kk++) {
                    new_Lambda[get_k_diag(k, kk, new_K)] = 
                         (kk == new_K - 1) ? R::rnorm(mu, sqrt(sigsq)) : Lambda[get_k_diag(k, kk, K)];
               }
          }
          
          // Compute the log prior term for Lambda values in the open table
          log_prior_term = 0.0;
          for (uword k = 0; k < new_K - 1; k++) {  // Only loop over interactions with the new cluster
               log_prior_term += R::dnorm(new_Lambda[get_k_diag(k, new_K - 1, new_K)], mu, sqrt(sigsq), 1);
          }
          
          for (uword j = 0; j < I; j++) {
               if (i == j) 
                    continue;
               
               y_index = (i < j) ? get_k(i, j, I) : get_k(j, i, I);
               lambda_index = get_k_diag(std::min(new_K - 1, Xi[j]), std::max(new_K - 1, Xi[j]), new_K);
               
               logprobs[new_K - 1] += R::dbinom(Y[y_index], 1, R::pnorm(new_Lambda[lambda_index], 0.0, 1.0, 1, 0), 1);
          }
          
          // Add the prior term to the log probability of the new cluster
          logprobs[K] += log_prior_term;
          
          // Sample new assignment
          uword new_assignment = wsample(exp(logprobs - logprobs.max()));
          
          // Temporarily remove the current assignment
          uword prev_cluster = Xi[i];
          cluster_sizes[prev_cluster]--;
          Xi[i] = new_assignment;
          
          // Update Xi and Lambda if a new cluster is created
          if (new_assignment == K) {  
               cluster_sizes.insert_rows(K, 1);
               Lambda = new_Lambda;
               K = new_K;
          }
          
          // Relabel Xi and reshape Lambda if a cluster becomes empty
          if (cluster_sizes[prev_cluster] == 0 && new_assignment != prev_cluster) {
               
               arma::uvec valid_clusters = find(cluster_sizes > 0);
               new_K = valid_clusters.n_elem;
               
               // Map old clusters to new indices
               arma::uvec old_to_new(K, fill::zeros);
               for (uword k = 0; k < new_K; k++) {
                    old_to_new[valid_clusters[k]] = k;
               }
               
               // Update Xi with the new cluster numbering
               for (uword j = 0; j < I; j++) {
                    Xi[j] = old_to_new[Xi[j]];
               }
               
               // Reshape Lambda to remove empty cluster rows and columns
               arma::vec new_Lambda(new_K * (new_K + 1) / 2, fill::zeros);
               
               for (uword k = 0; k < new_K; k++) {
                    for (uword kk = k; kk < new_K; kk++) {
                         k_old = valid_clusters[k];
                         kk_old = valid_clusters[kk];
                         new_Lambda[get_k_diag(k, kk, new_K)] = Lambda[get_k_diag(k_old, kk_old, K)];
                    }
               }
               
               Lambda = new_Lambda;
               K = new_K;
          }
     }
     
     return List::create(
          Named("Xi") = Xi,
          Named("Lambda") = Lambda,
          Named("K") = K
     );
}

// [[Rcpp::export]]
List sample_Lambda(const double& b, int n_tun_Lambda, double del2_Lambda, int n_Lambda, const int& n_burn, const double& mean_K, 
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
          double mix_rate = static_cast<double>(n_Lambda) / (b * 0.5 * mean_K * (mean_K + 1.0));
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
rowvec incidence_matrix0(const double& I, const double& B, const umat& Xi_chain) 
{
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
rowvec interaction_probs0(const double& I, const double& B, const List& Lambda_chain, const umat& Xi_chain) {
     const uword N = 0.5 * I * (I - 1.0);
     rowvec out(N, fill::zeros);
     const double scale_factor = 1.0 / B;
     
     urowvec Xi(I);
     rowvec Lambda;
     uword k, lambda_index;
     
     for (uword b = 0; b < B; b++) {
          
          // Get Lambda vector for sample b
          Lambda = as<rowvec>(Lambda_chain[b]);
          
          // Get Xi for sample b
          Xi = Xi_chain.row(b);
          
          // Now, determine K dynamically
          uword K = (-1.0 + std::sqrt(1.0 + 8.0 * Lambda.n_elem)) / 2.0;  // Solve 0.5 * K * (K + 1) = length(Lambda)
          
          for (uword i = 0; i < I - 1; i++) {
               for (uword ii = i + 1; ii < I; ii++) {
                    k = get_k(i, ii, I);
                    lambda_index = get_k_diag(std::min(Xi[i], Xi[ii]), std::max(Xi[i], Xi[ii]), K);
                    out(k) += R::dnorm(Lambda[lambda_index], 0.0, 1.0, 0) * scale_factor;
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
List WAIC(const double& I, const double& B, const vec& Y, const List& Lambda_chain, 
             const umat& Xi_chain) {
     uword m, lambda_index;
     double tmp, a_ib, a_ib_ssq, a_ib_sum;
     const double scale_factor = 1.0 / B;
     
     urowvec Xi(I);
     double lppd = 0.0, slp = 0.0, pWAIC2 = 0.0;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               m = get_k(i, ii, I);
               tmp = 0.0;
               a_ib_ssq = 0.0;
               a_ib_sum = 0.0;
               
               for (uword b = 0; b < B; b++) {
                    vec Lambda = Lambda_chain[b];  // Extract vector for iteration b
                    Xi = Xi_chain.row(b);
                    
                    uword K = arma::max(Xi) + 1;  // Infer K* from Xi (0-based indexing)
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
mat simulate_data(const double& I, const vec& Lambda, const uvec& Xi) {
     mat Y(I, I, fill::zeros);
     uword lambda_idx, K = arma::max(Xi) + 1;
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