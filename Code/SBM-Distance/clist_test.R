rm(list=ls(all=TRUE))

loc     <- "~/Dropbox/PAPERS/JNETS/"
dataset <- "test" 
model   <- "clist"

suppressMessages(suppressWarnings(library(Rcpp)))

sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
source   (paste0(loc, "code/rfunctions.R"))
source   (paste0(loc, "code/", model, "_simulate.R"))
source   (paste0(loc, "code/", model, "_functions.R"))

# GENERATE DATA
K      <- 5
Q      <- 2
I      <- K*30
zeta   <- 0
sigsq  <- 6
Xi     <- as.matrix(rep(1:K, rep(I/K, K)))

set.seed(1)
synda <- generate_data(I, K, Q, sigsq, zeta, Xi)

# load("C:/Users/Juan Camilo/Dropbox/PAPERS/JNETS/data/lazega_data.RData")

# MCMC
n_sams <- 100000
n_burn <- 100000
n_skip <- 1
ptm <- proc.time()
MCMC(Y, K, Q, n_sams, n_burn, n_skip, model, dataset, loc)
load(paste0(loc, "outs/", model, "_", dataset, "_", "samples.RData")) 
proc.time() - ptm

# CHAINS
plot_chain(THETA$loglik_chain, ll, main = "Log-likelihood")
plot_chain(THETA$Eta_chain, Eta, main = expression(eta), subset_indices = T)
plot_chain(THETA$zeta_chain, zeta, main = expression(zeta))
plot_chain(THETA$sigsq_chain, sigsq, main = expression(sigma^2))
# plot_ci(THETA$U_chain, U)
plot_chain0(THETA$omesq_chain, main = expression(omega^2))
plot_chain0(THETA$alpha_chain, main = expression(alpha))
plot_chain0(THETA$omega_chain, main = expression(omega))

# INCIDENCE MATRIX
incidence_matrix(THETA)

# INTERACTION PROBABILITIES
# interaction_probs(K, THETA, model, dataset, loc)
