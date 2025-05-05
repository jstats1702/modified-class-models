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
                    out(k) += R::dnorm(Lambda[lambda_index], 0.0, 1.0, 0) * scale_factor;
               }
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