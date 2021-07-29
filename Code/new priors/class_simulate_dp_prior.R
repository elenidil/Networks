
#rename the file to class_test_dp_prior
#install.packages(c("installr", "Rcpp"))
library(installr)
#install.Rtools()
rm(list=ls(all=TRUE))


library(Rcpp)
#evalCpp("1+1")

setwd("C:/Users/elena/Desktop/ESBM/code/")
loc     <- "C:/Users/elena/Desktop/ESBM/"
dataset <- "test_dp" 
model   <- "class"

suppressMessages(suppressWarnings(library(Rcpp)))
suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(MCMCpack)))
suppressMessages(suppressWarnings(library(mcmcplots)))

sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
sourceCpp(paste0(loc, "code/", model, "_functions_dp.cpp"))
sourceCpp(paste0(loc, "code/", model, "_simulate_dp_prior.cpp"))

#Note: the cpp code has to be written in cpp and compiled otherwise it doesn't seem to work in R



source(paste0(loc, "code/rfunctions.R")) #this works for all priors
#the class_simulate is used only to simulate the data - if I have already data is not needed
source(paste0(loc, "code/", model, "_simulate.R")) #those two need to be removed
source(paste0(loc, "code/", model, "_functions_dp.R")) #and copy paste the functions here
#or can make the file class_functions_dp.R which should be less messy
# LOAD DATA
K          <- 10
I          <- K*15
omega      <- rep(1,K)/K
mu_mu      <- -1
sig2_mu    <- 0.01
a_sig      <- 2 + 0.05^(-2)
b_sig      <- (a_sig-1)*1

set.seed(1234)
I <- 200
Ycube <- read.csv("C:/Users/elena/Desktop/ESBM/outs/sim10_new/Ycube.csv")
iP <- Y <- matrix(NA, I*(I-1)/2, 1)
for (i in 1:(I-1)) {
  for (ii in (i+1):I) {
    #iP[get_k(i, ii, I)] <- iPcube[i, ii]
    Y [get_k(i, ii, I)] <- Ycube [i, ii]
  }
}
image(as.matrix(Ycube))
#synda <- generate_data(I, K, omega, mu_mu, sig2_mu, a_sig, b_sig, model, dataset, loc)

#code to include more priors for the cluster assignments
#To choose the hyperparameter alpha_dp such that E[H] = 10
sigma_dp <- 0
H_dp <- Inf 
alpha_dp <- 2.55
V <- 100
expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = H_dp)

#those are probably not needed
urn_DP <- function(v_minus,alpha_PY){
  return(c(v_minus,alpha_PY))
}
if (prior=="DP"){
  urn<-function(v_minus){return(urn_DP(v_minus,alpha_PY))}
}
log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
#------
# MCMC
n_sams <- 2000
n_burn <- 2000
n_skip <- 1

MCMC(Y, K, n_sams, n_burn, n_skip, model, dataset, loc)
load(paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "samples.RData")) 

# CODA
plot_chain (THETA$Lambda_chain, Lambda, main = expression(lambda), subset_indices = T)
plot_chain (THETA$mu_chain, mu, main = expression(mu))
plot_chain (THETA$sigsq_chain, sigsq, main = expression(sigma^2))
plot_chain0(THETA$alpha_chain, main = expression(alpha))
plot_chain0(THETA$omega_chain, main = expression(omega))

# INCIDENCE MATRIX
incidence_matrix(THETA)

