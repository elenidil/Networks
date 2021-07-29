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


// This is a simple example of exporting a C++ function to R. You can source this function into 
//an R session using the Rcpp::sourceCpp function (or via the Source button on the editor toolbar).
//Learn more about Rcpp at:
//   http://www.rcpp.org/ //   http://adv-r.had.co.nz/Rcpp.html //   http://gallery.rcpp.org/

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.

uword get_k (const int& i, const int& ii, const double& I)
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
double loglik (const double& I, const double& K, const vec& Lambda, const uvec& Xi, const vec& Y)
{
  // Xi is a I x 1 vec whose indices are zero-based
  double out = 0.0;
  for (uword i = 0; i < I-1; i++) {
    for (uword ii = i+1; ii < I; ii++) {
      out += R::dbinom(Y[get_k(i, ii, I)], 1, expit(Lambda[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)]), 1);
    }
  }
  return(out);
}

// [[Rcpp::export]]
double sample_sigsq (const double& K, const double& a_sig, const double& b_sig, const double& mu, const vec& Lambda)
{
  return(1.0/R::rgamma(a_sig + K*(K+1.0)/4.0, 1.0/(b_sig + 0.5*accu(pow(Lambda - mu, 2)))));
}

// [[Rcpp::export]]
double sample_mu (const double& K, const double& mu_mu, const double& sig2_mu, const double& sigsq, const vec& Lambda)
{
  double v2 = 1.0/(1.0/sig2_mu + (0.5*K*(K+1.0))/sigsq);
  double m  = v2*(mu_mu/sig2_mu + accu(Lambda)/sigsq);
  return(R::rnorm(m, sqrt(v2)));
}


//#this used to be in cpp file where counting starts from 0, here if I let Xi as it is there are clusters 0
//#which causes errors when run in R and not in cpp, that's why I added +1 
//#Question: When it creates a new cluster so K <- K+1 it gives an error

// [[Rcpp::export]]

uvec sample_Xi_dp (const double& I, const double& K, const vec& Lambda, uvec Xi, const vec& Y, const vec& com_siz, const double& alpha_dp)
{
  // Xi contains uwords from 0 (class 1) to K-1 (class K)
  for (uword i = 0; i < I; i++) {
    vec logprobs(K, fill::zeros);
    for (uword k = 0; k < K; k++) {
      if (com_siz[k] > 0 )  logprobs[k] += log(com_siz[k]);
      if (com_siz[k] == 0 ) logprobs[k] += log(alpha_dp); 
      if (i < I-1) for (uword ii = i+1; ii < I; ii++) logprobs[k] += R::dbinom(Y[get_k(i, ii, I)], 1, expit(Lambda[get_k_diag(min(k, Xi[ii]), max(k, Xi[ii]), K)]), 1);
      if (i > 0  ) for (uword ii = 0;   ii < i; ii++) logprobs[k] += R::dbinom(Y[get_k(ii, i, I)], 1, expit(Lambda[get_k_diag(min(k, Xi[ii]), max(k, Xi[ii]), K)]), 1);
    }
    Xi[i] = wsample(exp(logprobs - max(logprobs)));
  }
  return(Xi);
}

// [[Rcpp::export]]
List sample_Lambda (const double& b, int n_tun_Lambda, double del2_Lambda, int n_Lambda, const int& n_burn, 
                    const double& I, const double& K, const double& sigsq, const double& mu, vec Lambda, const uvec& Xi, const vec& Y)
{
  
  // simplificar log verosimilitud
  uword m;
  double skl, nkl, lambda_c, lambda_p;
  for (uword k = 0; k < K; k++) {
    for (uword l = k; l < K; l++) {
      // suffitient statistics
      uvec ind_k = find(Xi == k);
      uvec ind_l = find(Xi == l);
      skl = 0.0;
      nkl = 0.0;
      for (uword i = 0; i < ind_k.n_elem; i++) {
        for (uword ii = 0; ii < ind_l.n_elem; ii++) {
          if (ind_k[i] < ind_l[ii]) {
            skl += (double)(Y[get_k(ind_k[i], ind_l[ii], I)]);
            nkl++;
          }
        }
      }
      // Metropolis step
      m = get_k_diag(k, l, K);
      lambda_c = Lambda[m];
      lambda_p = R::rnorm(lambda_c, sqrt(del2_Lambda));
      if (R::runif(0, 1) < exp(skl*(log(expit(lambda_p))-log(expit(lambda_c))) + (nkl-skl)*(log(1.0 - expit(lambda_p))-log(1.0 - expit(lambda_c))) + (-0.5/sigsq)*(pow(lambda_p - mu, 2)-pow(lambda_c - mu, 2)))) {
        Lambda[m] = lambda_p;
        n_Lambda++;
      }
    }
  }
  if (b < n_burn) {
    double mix_rate = n_Lambda/(b*0.5*K*(K+1.0));
    tunning(del2_Lambda, n_tun_Lambda, mix_rate, b);
  }
  return List::create(Named("Lambda")       = Lambda,
                      Named("del2_Lambda")  = del2_Lambda,
                      Named("n_Lambda")     = n_Lambda,
                      Named("n_tun_Lambda") = n_tun_Lambda);
}

// [[Rcpp::export]]
vec sample_Y (const double& I, const double& K, const vec& Lambda, const uvec& Xi, const vec& na_indices, vec Yna)
{
  // sample NA values in Y
  uword k;
  for (uword i = 0; i < I-1; i++) {
    for (uword ii = i+1; ii < I; ii++) {
      k = get_k(i, ii, I);
      if (na_indices[k] == true) Yna[k] = R::rbinom(1, expit(Lambda[get_k_diag(max(Xi[i], Xi[ii]), min(Xi[i], Xi[ii]), K)]));
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
rowvec interaction_probs0 (const double& I, const double& K, const double& B, const mat& Lambda_chain, const umat& Xi_chain) 
{
  urowvec Xi(I);
  rowvec Lambda(0.5*K*(K+1.0)), out(0.5*I*(I-1.0), fill::zeros);
  for (uword b = 0; b < B; b++) {
    Lambda = Lambda_chain.row(b);
    Xi = Xi_chain.row(b);
    for (uword i = 0; i < I-1; i++) {
      for (uword ii = i+1; ii < I; ii++) {
        out[get_k(i, ii, I)] += expit(Lambda[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)])/B;
      }
    }
  }
  return(out);
}

// [[Rcpp::export]]
List WAIC (const double& I, const double& K, const double& B, const vec& Y, const mat& Lambda_chain, const umat& Xi_chain)
{
  uword i, ii, m, b;
  double tmp, a_ib, a_ib_ssq, a_ib_sum;
  urowvec Xi(I);
  rowvec Lambda(0.5*K*(K+1.0));
  double lppd = 0.0, slp = 0.0, pWAIC2 = 0.0;
  for (i = 0; i < I-1; i++) {
    for (ii = i+1; ii < I; ii++) {
      m = get_k(i, ii, I);
      tmp = 0.0;
      a_ib_ssq = 0.0;
      a_ib_sum = 0.0;
      for (b = 0; b < B; b++) {
        Lambda = Lambda_chain.row(b);
        Xi = Xi_chain.row(b);
        a_ib = R::dbinom(Y[m], 1, expit(Lambda[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)]), 1);
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
mat simulate_data (const double& I, const double& K, const vec& Lambda, const uvec& Xi)
{
  mat Y(I, I, fill::zeros);
  for (uword i = 0; i < I-1; i++) {
    for (uword ii = i+1; ii < I; ii++) {
      Y.at(i, ii) = R::rbinom(1, expit(Lambda[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)]));
    }
  }
  Y = Y + Y.t();
  return(Y);
}
