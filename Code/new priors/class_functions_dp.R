get_hyperpars <- function()  #changed
{
  # hyper-parameters
  mu_mu   <- 0
  sig2_mu <- 3
  a_sig   <- 2
  b_sig   <- 3
  alpha_dp <- 2.55  #this has to be defined such that E[H]=10
  #a_alpha <- 1;  #b_alpha <- 1
  # MH parameters
  del2_Lambda  <- 1
  n_Lambda     <- 0
  n_tun_Lambda <- 100
  #del2_alpha   <- 1;  #n_alpha      <- 0;  #n_tun_alpha  <- 100
  list(mu_mu = mu_mu, sig2_mu = sig2_mu, a_sig = a_sig, b_sig = b_sig, alpha_dp = alpha_dp, 
       del2_Lambda = del2_Lambda, n_Lambda = n_Lambda, n_tun_Lambda = n_tun_Lambda)
}

get_initial_values <- function(I, K, hyps) #changed
{
  alpha_dp  <- 1 #SHOULD I CHANGE THAT?
  #omega  <- rep(alpha/K, K) #removed omega from list at the end
  # non random initial cluster assignment
  #Xi     <- as.matrix(c(rep(0, I%%K), rep(0:(K-1), rep(floor(I/K), K))))
  #random initial cluster assignment
  #It used to be 0:(K-1), but this causes problems when cpp is not used, so I'll change it to 1:K
  #maybe need to change it back if I use the cpp functions only or have 2 versions (one for R functions and one for cpp)
  Xi <- as.matrix(sample(0:(K-1), I, replace = T)) 
  sigsq  <- 1/rgamma(n = 1, shape = hyps$a_sig, rate = hyps$b_sig)
  mu     <- rnorm(n = 1, mean = hyps$mu_mu, sd = sqrt(hyps$sig2_mu))
  Lambda <- rnorm(n = K*(K+1)/2, mean = mu, sd = sqrt(sigsq))
  list(alpha_dp = alpha_dp, Xi = Xi, sigsq = sigsq, mu = mu, Lambda = Lambda)
}
get_chains_data <- function(I, K, n_sams)
{
  Lambda_chain <- matrix(NA, n_sams, K*(K+1)/2)
  Xi_chain     <- matrix(NA, n_sams, I)
  # mu_chain     <- matrix(NA, n_sams, 1);  # sigsq_chain  <- matrix(NA, n_sams, 1)
  # omega_chain  <- matrix(NA, n_sams, K);  # alpha_chain  <- matrix(NA, n_sams, 1)
  # loglik_chain <- matrix(NA, n_sams, 1)
  # list(Lambda_chain = Lambda_chain, mu_chain = mu_chain, sigsq_chain = sigsq_chain,
  #      Xi_chain = Xi_chain, omega_chain = omega_chain, alpha_chain = alpha_chain, 
  #      loglik_chain = loglik_chain)
  list(Lambda_chain = Lambda_chain, Xi_chain = Xi_chain)
}

#sample_Xi_dp <- function(I, K, Lambda, Xi, Y, com_siz, alpha_dp){
#  //#Xi_rver <- Xi + 1
#    for(i in 1:I){
#      logprobs <- rep(0, K)
#      for(k in 1: K){
#        if(is.na(com_siz[k])){
#          if(com_siz[k] > 0){ logprobs[k] = logprobs[k] + log(com_siz[k]) }
#          else if(com_siz[k] == 0){ logprobs[k] = logprobs[k] + log(alpha_dp)}
#        }
#        if(i < (I-1)) {for (ii in (i+2):I) { 
#          logprobs[k] = logprobs[k] + 
#            dbinom(Y[get_k(i,ii,I)], 1, expit(Lambda[get_k_diag(min(k, Xi[ii]+1), max(k, Xi[ii]+1), K)]), 1)} }
#        if(i > 1) {for (ii in 1:(i-1)) { 
#          logprobs[k] = logprobs[k] + 
#            dbinom(Y[get_k(ii,i,I)], 1, expit(Lambda[get_k_diag(min(k, Xi[ii]+1), max(k, Xi[ii]+1), K)]), 1)}  }
#      }
#      Xi[i] <- wsample(exp(logprobs - max(logprobs))) //#I should add -1 since I have the +1 above?
#    } 
#  Xi //#not sure if this is needed 
#} 
  

MCMC <- function(Y, K, n_sams, n_burn, n_skip, model, dataset, loc) #changed
{
  I       <- get_I(Y)
  THETA   <- get_chains_data(I, K, n_sams)
  ## hyper-parameters
  hyps    <- get_hyperpars()
  mu_mu   <- hyps$mu_mu
  sig2_mu <- hyps$sig2_mu
  a_sig   <- hyps$a_sig
  b_sig   <- hyps$b_sig
  alpha_dp <- hyps$alpha_dp
  #a_alpha <- hyps$a_alpha; b_alpha <- hyps$b_alpha
  ## MH parameters
  del2_Lambda  <- hyps$del2_Lambda
  n_Lambda     <- hyps$n_Lambda
  n_tun_Lambda <- hyps$n_tun_Lambda
  #del2_alpha   <- hyps$del2_alpha;  n_alpha<- hyps$n_alpha; n_tun_alpha  <- hyps$n_tun_alpha
  ## initial values
  tmp    <- get_initial_values(I, K, hyps)
  Lambda <- tmp$Lambda
  mu     <- tmp$mu
  sigsq  <- tmp$sigsq
  Xi     <- tmp$Xi
  #omega  <- tmp$omega
  #alpha_dp <- tmp$alpha_dp I think the one in hyper parameters is enough
  ## chains
  B <- n_burn + n_skip*n_sams
  n_disp <- floor(0.1*B)
  for (b in 1:B) 
  {
    tmp          <- sample_Lambda(b, n_tun_Lambda, del2_Lambda, n_Lambda, n_burn, I, K, sigsq, mu, Lambda, Xi, Y)
    Lambda       <- tmp$Lambda
    del2_Lambda  <- tmp$del2_Lambda
    n_Lambda     <- tmp$n_Lambda
    n_tun_Lambda <- tmp$n_tun_Lambda
    #I don't need to sample alpha_dp right now. Maybe later we'll change it (check again with ESBM paper)
    #tmp    <- sample_alpha(b, n_tun_alpha, del2_alpha, n_alpha, n_burn, K, a_alpha, b_alpha, alpha, omega)
    #alpha <- tmp$alpha; del2_alpha  <- tmp$del2_alpha; n_alpha <- tmp$n_alpha; n_tun_alpha <- tmp$n_tun_alpha
    mu           <- sample_mu(K, mu_mu, sig2_mu, sigsq, Lambda)
    sigsq        <- sample_sigsq(K, a_sig, b_sig, mu, Lambda)
    #create a vector with the sizes of each community (this can be done in Rcpp as well)
    com_siz <- rep(0,K)
    com_siz <- table(Xi)
    #added in the arguments of sample_Xi the com_siz and the alpha_dp
    Xi           <- sample_Xi_dp(I, K, Lambda, Xi, Y, com_siz, alpha_dp)
    #omega        <- sample_omega(K, alpha, Xi)
    # store
    if ((b > n_burn) & (b%%n_skip == 0)) {
      i <- (b - n_burn)/n_skip
      THETA$loglik_chain[i]  <- loglik(I, K, Lambda, Xi, Y)
      THETA$Xi_chain[i,]     <- c(Xi)
      THETA$Lambda_chain[i,] <- c(Lambda)
      # THETA$alpha_chain[i]   <- c(alpha); # THETA$mu_chain[i]      <- c(mu)
      # THETA$sigsq_chain[i]   <- c(sigsq); # THETA$omega_chain[i,]  <- c(omega)
    }
    # display        
    if (b%%n_disp == 0) {
      cat("Progress: ", dec(100*b/B, 1), "% complete", "\n",
          # "Lambda mixing rate = ", dec(100*n_Lambda/(b*0.5*K*(K+1)), 2), "%, del2_Lambda = ", dec(del2_Lambda, 5), "\n",
          # "alpha mixing  rate = ", dec(100*n_alpha/b,                2), "%, del2_alpha  = ", dec(del2_alpha,  5), "\n",
          "---------------------------------------------------", "\n", sep="")
    }
  }
  save(THETA, file = paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "samples.RData"))
}


YPPP <- function(Yna, na_indices, K, n_sams, n_burn, n_skip)  #changed
{
  I <- get_I(Yna)
  ## hyper-parameters
  hyps    <- get_hyperpars()
  mu_mu   <- hyps$mu_mu
  sig2_mu <- hyps$sig2_mu
  a_sig   <- hyps$a_sig
  b_sig   <- hyps$b_sig
  alpha_dp <- hyps$alpha_dp
  ## MH parameters
  del2_Lambda  <- hyps$del2_Lambda
  n_Lambda     <- hyps$n_Lambda
  n_tun_Lambda <- hyps$n_tun_Lambda
  ## initial values
  tmp    <- get_initial_values(I, K, hyps)
  Lambda <- tmp$Lambda
  mu     <- tmp$mu
  sigsq  <- tmp$sigsq
  Xi     <- tmp$Xi
  ## posterior predictive probabilities
  y_ppp <- rep(0, sum(na_indices))
  ## chains
  B <- n_burn + n_skip*n_sams
  for (b in 1:B) 
  {
    Yna          <- sample_Y(I, K, Lambda, Xi, na_indices, Yna)
    tmp          <- sample_Lambda(b, n_tun_Lambda, del2_Lambda, n_Lambda, n_burn, I, K, sigsq, mu, Lambda, Xi, Yna)
    Lambda       <- tmp$Lambda
    del2_Lambda  <- tmp$del2_Lambda
    n_Lambda     <- tmp$n_Lambda
    n_tun_Lambda <- tmp$n_tun_Lambda
    #deleted all the sample_alpha part
    mu           <- sample_mu(K, mu_mu, sig2_mu, sigsq, Lambda)
    sigsq        <- sample_sigsq(K, a_sig, b_sig, mu, Lambda)
    com_siz <- rep(0,K)
    com_siz <- table(Xi)
    Xi           <- sample_Xi_dp(I, K, Lambda, Xi, Y, com_siz, alpha_dp)
    #omega        <- sample_omega(K, alpha, Xi)
    # ppp
    if ((b > n_burn) & (b%%n_skip == 0)) y_ppp <- y_ppp + Yna[na_indices]/n_sams
  }
  y_ppp
}



GOF <- function(Y, THETA, model, dataset, loc)
{
  I <- get_I(Y)
  K <- (-1 + sqrt(1 + 8*ncol(THETA$Lambda_chain)))/2
  B <- nrow(THETA$Lambda_chain)
  gof <- WAIC(I, K, B, Y, THETA$Lambda_chain, THETA$Xi_chain)
  save(gof, file = paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "gof.RData"))
}

incidence_matrix <- function(THETA)
{
  I <- ncol(THETA$Xi_chain)
  B <- nrow(THETA$Xi_chain)
  K <- (-1 + sqrt(1 + 8*ncol(THETA$Lambda_chain)))/2
  A_vec <- incidence_matrix0(I, K, B, THETA$Xi_chain)
  A <- matrix(0, I, I)
  for (i in 1:(I-1)) for (ii in (i+1):I) A[i, ii] <- A_vec[get_k(i,ii,I)]
  A <- A + t(A)
  diag(A) <- 1
  ## plot
  pdf(paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "incidence_matrix.pdf"), height = 4, width = 4, pointsize = 10)
  par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
  heat.plot0(mat = A, show.grid = F, tick = 0, labs = NA, main = "Incidence matrix")
  dev.off()
}
# This is from simulate data folder, it is not needed if I already have the data, but if I want 
# to generate data to run simulations, I need to use it. The question is now, what do I need to use
# instead of omega to generate X. Everything else stays the same

generate_data <- function(I, K, mu_mu, sig2_mu, a_sig, b_sig, model, dataset, loc, plot_data = T)
{ #the function takes as argument omega, that is used to produced Xi
  suppressMessages(suppressWarnings(library(igraph)))
  ## parameters
  Xi     <- sort(sample(x = 1:K, size = I, replace = T)) #, prob = omega))
  mu     <- rnorm(n = 1, mean = mu_mu, sd = sqrt(sig2_mu))
  sigsq  <- 1/rgamma(n = 1, shape = a_sig, rate = b_sig)
  Lambda <- c(0.75, -1, -2,  -3, 0.75, -1, -2, 0.75, -1, 0.75)   # rnorm(n = K*(K+1)/2, mean = mu, sd = sqrt(sigsq))   
  ## Likelihood
  iPcube <- Ycube <- array(0, c(I, I))
  for(i in 1:(I-1)) {
    for(ii in (i+1):I) {
      iPcube[i, ii] <- iPcube[ii, i] <- expit(Lambda[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)])
      Ycube [i, ii] <- Ycube [ii, i] <- rbinom(n = 1, size = 1, prob = iPcube[i, ii])
    }
  }
  # as matrices
  iP <- Y <- matrix(NA, (I*(I-1)/2), 1)
  for (i in 1:(I-1)) {
    for (ii in (i+1):I) {
      iP[get_k(i, ii, I)] <- iPcube[i, ii]
      Y [get_k(i, ii, I)] <- Ycube [i, ii]
    }
  }
  ## plots
  if (plot_data) {
    pdf(paste0(loc, "outs/", dataset, "/", model, "_", dataset,"_adjacency.pdf"), height = 1*4, width = 1*4, pointsize = 10)
    par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
    adjacency.plot(mat = Ycube,  labs = NA, show.grid = F, tick = 0, main = "Adjacency Matrix")
    dev.off()
    
    pdf(paste0(loc, "outs/", dataset, "/", model, "_", dataset,"_probabilities.pdf"), height = 1*4, width = 1*4, pointsize = 10)
    par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
    heat.plot0    (mat = iPcube, labs = NA, show.grid = F, tick = 0, main = "Interaction Probabilities")
    dev.off()
    
    pdf(paste0(loc, "outs/", dataset, "/", model, "_", dataset,"_graph.pdf"), height = 1*4, width = 1*4, pointsize = 10)
    par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
    g <- graph.adjacency(adjmatrix = Ycube, mode = "undirected", diag = F)
    plot(g, vertex.color = "lightblue", edge.color = "gray50", vertex.label = NA, vertex.size = 3)
    title("Graph")
    box()
    dev.off()
  }
  ## load
  Xi     <-  Xi-1
  Xi     <<- Xi
  #omega  <<- omega
  mu     <<- mu
  sigsq  <<- sigsq
  Lambda <<- Lambda
  Y      <<- Y
  Ycube  <<- Ycube
  ll     <<- loglik(I, K, Lambda, Xi, Y)
  ## return
  list(Y = Y, Ycube = Ycube, Xi = Xi, mu = mu, sigsq = sigsq, Lambda = Lambda, ll = ll) #omega = omega,
}



GOF <- function(Y, THETA, model, dataset, loc)
{
  I <- get_I(Y)
  K <- (-1 + sqrt(1 + 8*ncol(THETA$Lambda_chain)))/2
  B <- nrow(THETA$Lambda_chain)
  gof <- WAIC(I, K, B, Y, THETA$Lambda_chain, THETA$Xi_chain)
  save(gof, file = paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "gof.RData"))
}

incidence_matrix <- function(THETA)
{
  I <- ncol(THETA$Xi_chain)
  B <- nrow(THETA$Xi_chain)
  K <- (-1 + sqrt(1 + 8*ncol(THETA$Lambda_chain)))/2
  A_vec <- incidence_matrix0(I, K, B, THETA$Xi_chain)
  A <- matrix(0, I, I)
  for (i in 1:(I-1)) for (ii in (i+1):I) A[i, ii] <- A_vec[get_k(i,ii,I)]
  A <- A + t(A)
  diag(A) <- 1
  ## plot
  pdf(paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "incidence_matrix.pdf"), height = 4, width = 4, pointsize = 10)
  par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
  heat.plot0(mat = A, show.grid = F, tick = 0, labs = NA, main = "Incidence matrix")
  dev.off()
}