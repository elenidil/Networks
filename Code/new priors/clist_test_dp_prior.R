#rename the file to clist_test_dp_prior
#install.packages(c("installr", "Rcpp"))
library(installr)
#install.Rtools()
rm(list=ls(all=TRUE))


library(Rcpp)

setwd("C:/Users/elena/Desktop/ESBM/code/")
loc     <- "C:/Users/elena/Desktop/ESBM/"
dataset <- "test_dp" 
model   <- "clist"

suppressMessages(suppressWarnings(library(Rcpp)))
suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(MCMCpack)))
suppressMessages(suppressWarnings(library(mcmcplots)))

#sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
sourceCpp(paste0(loc, "code/", model, "_functions_dp1.cpp"))
#sourceCpp(paste0(loc, "code/", model, "_functions_dp_prior.cpp"))

source(paste0(loc, "code/rfunctions.R")) #this works for all priors
#the class_simulate is used only to simulate the data - if I have already data is not needed
source(paste0(loc, "code/", model, "_simulate.R")) #those two need to be removed
source(paste0(loc, "code/", model, "_functions_dp.R")) #and copy paste the functions here
#or can make the file class_functions_dp.R which should be less messy
# LOAD/GENERATE DATA
#K          <- 5; I          <- K*15
#omega      <- rep(1,K)/K; mu_mu      <- -1; sig2_mu    <- 0.01
#a_sig      <- 2 + 0.05^(-2); b_sig      <- (a_sig-1)*1
#synda <- generate_data(I, K, omega, mu_mu, sig2_mu, a_sig, b_sig, model, dataset, loc)

set.seed(1234)
I <- 200
I <- 80
Ycube <- read.csv("C:/Users/elena/Desktop/ESBM/outs/sim10_new/Ycube.csv")
Ycube <- read.csv("C:/Users/elena/Desktop/ESBM/code/simulated data sets/network_1_80.csv")
Ycube <- Ycube[,-1]
iP <- Y <- matrix(NA, I*(I-1)/2, 1)
for (i in 1:(I-1)) {
  for (ii in (i+1):I) {
    #iP[get_k(i, ii, I)] <- iPcube[i, ii]
    Y [get_k(i, ii, I)] <- Ycube [i, ii]
  }
}
image(as.matrix(Ycube))

#code to include more priors for the cluster assignments
#To choose the hyperparameter alpha_dp such that E[H] = 10
sigma_dp <- 0
H_dp <- Inf 
alpha_dp <- 2.55  #this has to be changed in function get_hyperpars on file _fucntions_dp.R
#when K=20 alpha_dp is 7.2, when K=10 alpha_dp = 2.55, when K=5 alpha_dp = 0.945
V <- 100
expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = H_dp)

#------
library(tictoc)
# MCMC
n_sams <- 25000
n_burn <- 25000
n_skip <- 1
K <- 20
Q <- 2
tic()
MCMC(Y, K, Q, n_sams, n_burn, n_skip, model, dataset, loc)
toc()
#436 sec for 80x80 25k iter K=10 (true is 5)
load(paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "samples.RData")) 

# CODA
plot_chain (THETA$Lambda_chain, Lambda, main = expression(lambda), subset_indices = T)
plot_chain (THETA$mu_chain, mu, main = expression(mu))
plot_chain (THETA$sigsq_chain, sigsq, main = expression(sigma^2))
plot_chain0(THETA$alpha_chain, main = expression(alpha))
plot_chain0(THETA$omega_chain, main = expression(omega))

# INCIDENCE MATRIX
incidence_matrix(THETA)

#Note: For the 80x80 with 5 communities when I choose K = 5 and set the alpha_dp parameter such that
#E[H]=5 then the code doesn't predict anything. For K = 10 and alpha_dp=2.55 it works fine
