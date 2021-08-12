//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

uword get_k (const int& i, const int& ii, const double& I)
{
  return(0.5*(I*(I-1.0) - (I-i)*(I-i-1.0)) + ii - i - 1.0);
} 

uword get_k_test (const int i, const int ii, const double I)
{
  return(0.5*(I*(I-1.0) - (I-i)*(I-i-1.0)) + ii - i - 1.0);
} 
uword get_k_diag (const int& k, const int& kk, const double& K)
{
  return(K*k + kk - 0.5*k*(k+1.0));
}

double expit (const double& x) 
{
  return(1.0/(1.0 + exp(-x)));
}

void tunning (double& del2, int& n_tun, const double& mix_rate, const int& b)
{ 
  // tunning paramter calibration
  double mr = 0.35, eps = 0.05, diff = mix_rate - mr;
  if ( b%n_tun == 0 ) {
    if ( abs(diff) > eps ) {
      double tmp = del2, cont = 1.0;
      do {
        tmp = del2 + (0.1/cont) * diff;
        cont++;
      } while ( (tmp <= 0.0) && (cont <= 1000000.0) );
      if ( tmp > 0.0 ) del2 = tmp;
      n_tun = 100;
    } else {
      n_tun += 100;
    }
  }
}

uword wsample (vec probs) 
{
  // weighted sampling according to probs
  // return indices fron 0 to n - 1, where n = probs.n_elem
  probs = probs/accu(probs);
  double u = R::runif(0.0, 1.0);
  if (u < probs(0)) {
    return 0;
  } else {
    uword n = probs.n_elem, i, out = 0;
    vec probsum = cumsum(probs);
    for (i = 1; i < n; i++) {  
      if ((probsum(i-1) <= u) && (u < probsum(i))) {
        out = i;
        goto endloop;
      }
    }
    endloop:
      return( out );
  }
}
vec rdirichlet (const vec& alpha) 
{
  uword K = alpha.n_elem;
  vec out(K);
  for (uword k = 0; k < K; k++) out[k] = R::rgamma(alpha[k], 1.0);
  return(out/accu(out));
}

// [[Rcpp::export]]
double loglik (const double& I, const double& K, const vec& Eta, const uvec& Xi, const vec& Y)
{
  // Xi is a I x 1 vec whose indices are zero-based
  double out = 0.0;
  for (uword i = 0; i < I-1; i++) {
    for (uword ii = i+1; ii < I; ii++) {
      out += R::dbinom(Y[get_k(i, ii, I)], 1, expit(Eta[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)]), 1);
    }
  }
  return(out);
}

// [[Rcpp::export]]
vec get_Eta (const double& K, const double& zeta, const mat& U)
{
  // compute eta_{k,\ell}
  vec out(0.5*K*(K+1.0));
  for (uword k = 0; k < K; k++) {
    for (uword l = k; l < K; l++) {
      out[get_k_diag(k, l, K)] = zeta - norm(U.row(k) - U.row(l));
    }
  }
  return(out);
}

double lfcd_U (const rowvec& x, const uword& k, const double& I, const double& K, const double& sigsq, const double& zeta, mat U, const uvec& Xi, const vec& Y)
{
  double out = -pow(norm(x), 2)/(2.0*sigsq);
  U.row(k) = x;
  vec Eta = get_Eta(K, zeta, U);
  for (uword i = 0; i < I-1; i++) {
    for (uword ii = i+1; ii < I; ii++) {
      if ((Xi[i] == k) || (Xi[ii] == k)) {
        out += R::dbinom(Y[get_k(i, ii, I)], 1, expit(Eta[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)]), 1);
      }
    }
  }
  return(out);
}

// [[Rcpp::export]]
List sample_U (const double& b, int n_tun_U, double del2_U, int n_U, const int& n_burn, 
               const double& I, const double& K, const double& Q, const double& sigsq, const double& zeta, mat U, const uvec& Xi, const vec& Y)
{
  // Metropolis step
  rowvec u_c(Q), u_p(Q);
  for (uword k = 0; k < K; k++) {
    u_c = U.row(k);
    u_p = u_c + sqrt(del2_U)*randn<rowvec>(Q);
    if (R::runif(0, 1) < exp(lfcd_U(u_p, k, I, K, sigsq, zeta, U, Xi, Y) - lfcd_U(u_c, k, I, K, sigsq, zeta, U, Xi, Y))) {
      U.row(k) = u_p;
      n_U++;
    }
  }
  if (b < n_burn ) {
    double mix_rate = n_U/(b*K);
    tunning(del2_U, n_tun_U, mix_rate, b);
  }
  return List::create(Named("U")       = U,
                      Named("del2_U")  = del2_U,
                      Named("n_U")     = n_U,
                      Named("n_tun_U") = n_tun_U);
}

// [[Rcpp::export]]
List sample_zeta (const double& b, int n_tun_zeta, double del2_zeta, int n_zeta, const int& n_burn, 
                  const double& I, const double& K, const double& omesq, double zeta, const mat& U, const uvec& Xi, const vec& Y)
{
  // Metropolis step
  double zeta_p = R::rnorm(zeta, sqrt(del2_zeta));
  if (R::runif(0, 1) < exp(loglik(I, K, get_Eta(K, zeta_p, U), Xi, Y) - loglik(I, K, get_Eta(K, zeta, U), Xi, Y) - (pow(zeta_p, 2) - pow(zeta, 2))/(2.0*omesq))) {
    zeta = zeta_p;
    n_zeta++;
  }
  if (b < n_burn) {
    double mix_rate = n_zeta/b;
    tunning(del2_zeta, n_tun_zeta, mix_rate, b);
  }
  return List::create(Named("zeta")       = zeta,
                      Named("del2_zeta")  = del2_zeta,
                      Named("n_zeta")     = n_zeta,
                      Named("n_tun_zeta") = n_tun_zeta);
}

// [[Rcpp::export]]   
uvec sample_Xi_dp (const double& I, const double& K, const vec& Eta, uvec Xi, const vec& Y, const vec& com_siz, const double& alpha_dp)
{
  // Xi contains uwords from 0 (class 1) to K-1 (class K)
  for (uword i = 0; i < I; i++) {
    vec logprobs(K, fill::zeros);
    for (uword k = 0; k < K; k++) {
      if (com_siz[k] > 0 )  logprobs[k] += log(com_siz[k]);
      if (com_siz[k] == 0 ) logprobs[k] += log(alpha_dp);
      if (i < I-1) for (uword ii = i+1; ii < I; ii++) logprobs[k] += R::dbinom(Y[get_k(i, ii, I)], 1, expit(Eta[get_k_diag(min(k, Xi[ii]), max(k, Xi[ii]), K)]), 1);
      if (i > 0  ) for (uword ii = 0;   ii < i; ii++) logprobs[k] += R::dbinom(Y[get_k(ii, i, I)], 1, expit(Eta[get_k_diag(min(k, Xi[ii]), max(k, Xi[ii]), K)]), 1);
    }
    Xi[i] = wsample(exp(logprobs - max(logprobs)));
  }
  return(Xi);
}

// [[Rcpp::export]]  sample_omega THIS IS NOT NEEDED ALONG WITH sample_alpha
// maybe I can use later the sample_alpha when we want it to be sampled instead of constant

// [[Rcpp::export]]
double sample_sigsq (const double& K, const double &Q, const double& a_sig, const double& b_sig, const mat& U)
{
  return(1.0/R::rgamma(a_sig + 0.5*K*Q, 1.0/(b_sig + 0.5*accu(pow(U, 2)))));
}

// [[Rcpp::export]]
double sample_omesq (const double& a_ome, const double& b_ome, const double& zeta)
{
  return(1.0/R::rgamma(a_ome + 0.5, 1.0/(b_ome + 0.5*pow(zeta, 2))));
}


// [[Rcpp::export]]
vec sample_Y (const double& I, const double& K, const vec& Eta, const uvec& Xi, const vec& na_indices, vec Yna)
{
  // sample NA values in Y
  uword k;
  for (uword i = 0; i < I-1; i++) {
    for (uword ii = i+1; ii < I; ii++) {
      k = get_k(i, ii, I);
      if (na_indices[k] == true) Yna[k] = R::rbinom(1, expit(Eta[get_k_diag(max(Xi[i], Xi[ii]), min(Xi[i], Xi[ii]), K)]));
    }
  }
  return(Yna);
}
// [[Rcpp::export]]
rowvec incidence_matrix0 (const double& I, const double& K, const double& B, const umat& Xi_chain) 
{
  rowvec out(0.5*I*(I-1.0), fill::zeros);
  for (uword b = 0; b < B; b++) {
    for (uword i = 0; i < I-1; i++) {
      for (uword ii = i+1; ii < I; ii++) {
        if ( Xi_chain.at(b, i) == Xi_chain.at(b, ii) ) out[get_k(i, ii, I)] += 1.0/B;
      }
    }
  }
  return(out);
} 

// [[Rcpp::export]]
rowvec interaction_probs0 (const double& I, const double& K, const double& B, const mat& Eta_chain, const umat& Xi_chain) 
{
  urowvec Xi(I);
  rowvec Eta(0.5*K*(K+1.0)), out(I*(I-1.0)/2.0, fill::zeros);
  for (uword b = 0; b < B; b++) {
    Eta = Eta_chain.row(b);
    Xi  = Xi_chain.row(b);
    for (uword i = 0; i < I-1; i++) {
      for (uword ii = i+1; ii < I; ii++) {
        out[get_k(i, ii, I)] += expit(Eta[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)])/B;
      }
    }
  }
  return(out);
}

// [[Rcpp::export]]
List WAIC (const double& I, const double& K, const double& B, const vec& Y, const mat& Eta_chain, const umat& Xi_chain)
{
  uword i, ii, m, b;
  double tmp, a_ib, a_ib_ssq, a_ib_sum;{}
  urowvec Xi(I);
  rowvec Eta(0.5*K*(K+1.0));
  double lppd = 0.0, slp = 0.0, pWAIC2 = 0.0;
  for (i = 0; i < I-1; i++) {
    for (ii = i+1; ii < I; ii++) {
      m = get_k(i, ii, I);
      tmp = 0.0;
      a_ib_ssq = 0.0;
      a_ib_sum = 0.0;
      for (b = 0; b < B; b++) {
        Eta = Eta_chain.row(b);
        Xi = Xi_chain.row(b);
        a_ib = R::dbinom(Y[m], 1, expit(Eta[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)]), 1);
        // waic 1
        tmp += exp(a_ib)/B;
        slp += a_ib/B;
        // waic 2
        a_ib_ssq += pow(a_ib, 2);
        a_ib_sum += a_ib;
      }
      lppd   += log(tmp);
      pWAIC2 += (a_ib_ssq - B*pow(a_ib_sum/B, 2))/(B - 1.0);
    }
  }
  double pWAIC1 =  2.0*lppd - 2.0*slp;
  double waic1  = -2.0*lppd + 2.0*pWAIC1;
  double waic2  = -2.0*lppd + 2.0*pWAIC2;
  return List::create(Named("lppd")   = lppd,
                      Named("pWAIC1") = pWAIC1,
                      Named("pWAIC2") = pWAIC2,
                      Named("waic1")  = waic1,
                      Named("waic2")  = waic2);
}

// [[Rcpp::export]]
mat simulate_data (const double& I, const double& K, const vec& Eta, const uvec& Xi)
{
  mat Y(I, I, fill::zeros);
  for (uword i = 0; i < I-1; i++) {
    for (uword ii = i+1; ii < I; ii++) {
      Y.at(i, ii) = R::rbinom(1, expit(Eta[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)]));
    }
  }
  Y = Y + Y.t();
  return(Y);
}

