rm(list=ls(all=TRUE))

loc     <- "C:/Users/Juan Camilo/Dropbox/PAPERS/NETS/"
dataset <- "test" 
model   <- "class"

suppressMessages(suppressWarnings(library(Rcpp)))
suppressMessages(suppressWarnings(library(mvtnorm)))
suppressMessages(suppressWarnings(library(MCMCpack)))
suppressMessages(suppressWarnings(library(mcmcplots)))

sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))

source(paste0(loc, "code/rfunctions.R"))
source(paste0(loc, "code/", model, "_simulate.R"))
source(paste0(loc, "code/", model, "_functions.R"))

# LOAD DATA
K          <- 4
I          <- K*15
omega      <- rep(1,K)/K
mu_mu      <- -1
sig2_mu    <- 0.01
a_sig      <- 2 + 0.05^(-2)
b_sig      <- (a_sig-1)*1

set.seed(1234)
synda <- generate_data(I, K, omega, mu_mu, sig2_mu, a_sig, b_sig, model, dataset, loc)

# MCMC
n_sams <- 25000
n_burn <- 25000
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

# INTERACTION PROBABILITIES
# interaction_probs(K, THETA, model, dataset, loc)