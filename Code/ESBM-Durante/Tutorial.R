setwd('~/Documents/ESBM-master/Data and Codes')

rm(list=ls())
source("esbm.R")
library(Rcpp)
Rcpp::sourceCpp('stirling.cpp')

# This is my data with 20 communities of size 5
load("network_1O.RData")
V <- dim(Y)[1]

# Note that Y must have diagonal equal to 0
diag(Y)

# save for C++ code
write(Y, file="Ytest100_1O.txt", ncolumns=1)

rep <- "1O"
pdf(sprintf("adj%s_%s.pdf", V, rep)) 
par(mfrow=c(1,1),mar=c(2,2,1,0.5))
image(seq(1,V),seq(1,V),Y,col=c("white",rev(heat.colors(100))),axes=F)	
box()
dev.off()

# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------
sigma_dp <- 0   
H_dp <- Inf 
alpha_dp <- 15#2.55#15
expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = H_dp)

# ------------------------------------
# PITMAN-YOR PROCESS
# ------------------------------------
sigma_py <- 0.575
H_py <- Inf 
alpha_py <- 1.6
expected_cl_py(V, sigma = sigma_py, theta = alpha_py, H = H_py)

######## Run MCMC ########################## 

N_iter <- 20000
V <- dim(Y)[1]
my_seed <- 1
my_z <- c(1:V)

V <- 100
Kinit <- 10
set.seed(1234)
xiinit <- sample(1:Kinit, V, replace=T)
my_z <- as.numeric(factor(xiinit, labels = seq(1,length(unique(xiinit)))))

# save for C++ code 
write(as.matrix(my_z), file="xiinit100_1O.txt", ncolumns=1)


# ------------------------------------
# DIRICHLET PROCESS (CRP)
# ------------------------------------

my_prior <- "DP"
Z_DP <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_dp, sigma_dp)

burn <- 1:15000
z <- Z_DP[,-burn]
table(z[,5000])

# quantiles for number of non-empty clusters
burn_in <- 5000
summary(apply(Z_DP[,(burn_in+1):N_iter],2,max))
quantile(apply(Z_DP[,(burn_in+1):N_iter],2,max))[c(2:4)]

pairw <- pr_cc(z)

pdf(sprintf("pairw%s_%s_ranK10.pdf", V, rep))
par(mfrow=c(1,1),mar=c(2,2,2,0.5))
image(seq(1,V),seq(1,V), pairw,col=c("white",rev(heat.colors(100))),axes=F)
box()	
dev.off()

pp   <- pairw
start <- ceiling(order(lower.tri(pp,diag=F)*pp,decreasing=T)[ 1 ]/V) 
used   <- reorder.pairwise.alt(pp,start,V) 

pdf(sprintf("pairw%s_%s_ordered_ranK10.pdf", V, rep)) 
par(mfrow=c(1,1),mar=c(2,2,1,0.5))
image(seq(1,V),seq(1,V),pairw[used, used],col=c("white",rev(heat.colors(100))),axes=F)	
box()
dev.off()

#### r and p for ESC ########

lmus <- function(s, r, p) {
	lgm <- r*log(1-p)-log(1-(1-p)^r)
	lgm + lgamma(s + r) + s * log(p) - lgamma(s) - lgamma(s + 1)
}

s <- 1:10
r <- 50
p <- 0.5
(p*r)/(1-p)

mus <- exp(lmus(s, r, p))
plot(s, mus)
