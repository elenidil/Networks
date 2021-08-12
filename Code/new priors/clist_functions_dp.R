get_hyperpars <- function(K, Q)  #changed
{
  a_ome <- 2
  b_ome <- (a_ome-1)*3
  a_sig <- 2
  b_sig <- (a_sig-1)*(pi^(Q/2)/exp(lgamma(Q/2+1))*K^(2/Q))
  alpha_dp <- 7.2
  #a_alpha <- 1; b_alpha <- 1
  # MH parameters
  del2_U      <- 0.1
  n_U         <- 0
  n_tun_U     <- 100
  del2_zeta   <- 0.1
  n_zeta      <- 0
  n_tun_zeta  <- 100
  #del2_alpha  <- 1; n_alpha     <- 0; n_tun_alpha <- 100
  # return
  list(a_sig = a_sig, b_sig = b_sig, a_ome = a_ome, b_ome = b_ome, alpha_dp = alpha_dp, del2_U = del2_U, n_U = n_U, 
       n_tun_U = n_tun_U, del2_zeta = del2_zeta, n_zeta = n_zeta, n_tun_zeta = n_tun_zeta)
      # del2_alpha = del2_alpha, n_alpha = n_alpha, n_tun_alpha = n_tun_alpha, a_alpha = a_alpha, b_alpha = b_alpha)
}

get_initial_values <- function(I, K, Q, hyps) #changed
{
  #alpha_dp <- 1  is set in hyperparameters
  #omega <- rep(alpha/K, K)
  Xi    <- as.matrix(c(rep(0, I%%K), rep(0:(K-1), rep(floor(I/K), K))))
  #Xi <- as.matrix(sample(0:(K-1), I, replace = T)) 
  omesq <- 1/rgamma(n = 1, shape = hyps$a_ome, rate = hyps$b_ome)
  sigsq <- 1/rgamma(n = 1, shape = hyps$a_sig, rate = hyps$b_sig)
  zeta  <- rnorm(n = 1, mean = 0, sd = sqrt(omesq))
  U     <- matrix(rnorm(n = K*Q, mean = 0, sd = sqrt(sigsq)), K, Q)
  list(Xi = Xi, omesq = omesq, sigsq = sigsq, zeta = zeta, U = U) # omega = omega,alpha = alpha,
}

get_chains_data <- function(I, K, Q, n_sams)
{
  Eta_chain    <- matrix(NA, n_sams, 0.5*K*(K+1))
  Xi_chain     <- matrix(NA, n_sams, I)
  # U_chain      <- matrix(NA, n_sams, K*Q)
  # zeta_chain   <- matrix(NA, n_sams, 1)
  # sigsq_chain  <- matrix(NA, n_sams, 1)
  # omesq_chain  <- matrix(NA, n_sams, 1)
  # omega_chain  <- matrix(NA, n_sams, K)
  # alpha_chain  <- matrix(NA, n_sams, 1)
  loglik_chain <- matrix(NA, n_sams, 1)
  # list(alpha_chain = alpha_chain, omega_chain = omega_chain, Xi_chain = Xi_chain, 
  #    omesq_chain = omesq_chain, sigsq_chain = sigsq_chain, zeta_chain = zeta_chain, 
  #     U_chain = U_chain, Eta_chain = Eta_chain, loglik_chain = loglik_chain)
  list(Xi_chain = Xi_chain, Eta_chain = Eta_chain, loglik_chain = loglik_chain)
}

MCMC <- function(Y, K, Q, n_sams, n_burn, n_skip, model, dataset, loc) #changed
{
  I       <- get_I(Y) #(1 + sqrt(1 + 8*nrow(Y)))/2
  THETA   <- get_chains_data(I, K, Q, n_sams) #initializes Eta_chain and 
  #Xi_chain as two matrices of size n_sams* K(K+1)/2 for Eta and n_sams*I for Xi
  ## hyper-parameters
  hyps    <- get_hyperpars(K, Q)
  #a_alpha <- hyps$a_alpha;a_alpha <- hyps$b_alpha
  a_ome   <- hyps$a_ome;  b_ome   <- hyps$b_ome
  a_sig   <- hyps$a_sig;  b_sig   <- hyps$b_sig
  alpha_dp <- hyps$alpha_dp
  ## MH parameters
  del2_U      <- hyps$del2_U
  n_U         <- hyps$n_U
  n_tun_U     <- hyps$n_tun_U
  del2_zeta   <- hyps$del2_zeta
  n_zeta      <- hyps$n_zeta
  n_tun_zeta  <- hyps$n_tun_zeta
  #del2_alpha  <- hyps$del2_alpha;  n_alpha     <- hyps$n_alpha;  n_tun_alpha <- hyps$n_tun_alpha
  ## initial values
  tmp   <- get_initial_values(I, K, Q, hyps)
  U     <- tmp$U
  zeta  <- tmp$zeta
  Xi    <- tmp$Xi
  #omega <- tmp$omega;  alpha <- tmp$alpha
  sigsq <- tmp$sigsq
  omesq <- tmp$omesq
  ## chains
  B <- n_burn + n_skip*n_sams
  n_disp <- floor(0.01*B)
  for (b in 1:B) 
  {
    Eta         <- get_Eta(K, zeta, U)
    com_siz <- rep(0,K)
    com_siz <- table(Xi)
    Xi          <- sample_Xi_dp(I, K, Eta, Xi, Y, com_siz, alpha_dp)
    #omega       <- sample_omega(K, alpha, Xi)
    sigsq       <- sample_sigsq(K, Q, a_sig, b_sig, U)
    omesq       <- sample_omesq(a_ome, b_ome, zeta)
    #Sampling alpha is not needed right now
    #tmp         <- sample_alpha(b, n_tun_alpha, del2_alpha, n_alpha, n_burn, K, a_alpha, b_alpha, alpha, omega)
    #alpha       <- tmp$alpha; del2_alpha  <- tmp$del2_alpha; n_alpha     <- tmp$n_alpha; n_tun_alpha <- tmp$n_tun_alpha
    tmp         <- sample_U(b, n_tun_U, del2_U, n_U, n_burn, I, K, Q, sigsq, zeta, U, Xi, Y)
    U           <- tmp$U
    del2_U      <- tmp$del2_U
    n_U         <- tmp$n_U
    n_tun_U     <- tmp$n_tun_U
    tmp         <- sample_zeta(b, n_tun_zeta, del2_zeta, n_zeta, n_burn, I, K, omesq, zeta, U, Xi, Y)
    zeta        <- tmp$zeta
    del2_zeta   <- tmp$del2_zeta
    n_zeta      <- tmp$n_zeta
    n_tun_zeta  <- tmp$n_tun_zeta
    # store
    if ((b > n_burn) & (b%%n_skip == 0)) {
      i <- (b - n_burn)/n_skip
      THETA$loglik_chain[i] <- loglik(I, K, Eta, Xi, Y)
      THETA$Xi_chain[i,]    <- c(Xi)
      THETA$Eta_chain[i,]   <- c(Eta)
      # THETA$alpha_chain[i]  <- c(alpha);      # THETA$U_chain[i,]     <- c(U)
      # THETA$zeta_chain[i]   <- c(zeta);      # THETA$sigsq_chain[i]  <- c(sigsq)
      # THETA$omesq_chain[i]  <- c(omesq);      # THETA$omega_chain[i,] <- c(omega)
    }
    # display        
    # if (b%%n_disp == 0) {
    #         cat("Progress: ", dec(100*b/B, 1), "% complete", "\n",
    #             "U     mixing rate = ", dec(100*n_U/(b*K), 2), "%, del2_U     = ", dec(del2_U,     5),  "\n",
    #             "zeta  mixing rate = ", dec(100*n_zeta/b,  2), "%, del2_zeta  = ", dec(del2_zeta,  5),  "\n",
    #             "alpha mixing rate = ", dec(100*n_alpha/b, 2), "%, del2_alpha = ", dec(del2_alpha, 5),  "\n",
    #             "\n", sep="")
    # }
  }
  save(THETA, file = paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "samples.RData"))  
}

YPPP <- function(Yna, na_indices, K, Q, n_sams, n_burn, n_skip)  #changed
{
  I       <- get_I(Yna)
  ## hyper-parameters
  hyps    <- get_hyperpars(K, Q)
  alpha_dp <- hyps$alpha_dp
  a_ome   <- hyps$a_ome
  b_ome   <- hyps$b_ome
  a_sig   <- hyps$a_sig
  b_sig   <- hyps$b_sig
  ## MH parameters
  del2_U      <- hyps$del2_U
  n_U         <- hyps$n_U
  n_tun_U     <- hyps$n_tun_U
  del2_zeta   <- hyps$del2_zeta
  n_zeta      <- hyps$n_zeta
  n_tun_zeta  <- hyps$n_tun_zeta
  #del2_alpha  <- hyps$del2_alpha;  n_alpha     <- hyps$n_alpha;  n_tun_alpha <- hyps$n_tun_alpha
  ## initial values
  tmp   <- get_initial_values(I, K, Q, hyps)
  U     <- tmp$U
  zeta  <- tmp$zeta
  Xi    <- tmp$Xi
  #omega <- tmp$omega; alpha <- tmp$alpha
  sigsq <- tmp$sigsq
  omesq <- tmp$omesq
  ## posterior predictive probabilities
  y_ppp <- rep(0, sum(na_indices))
  ## chains
  B <- n_burn + n_skip*n_sams
  for (b in 1:B) 
  {
    Eta         <- get_Eta(K, zeta, U)
    Yna         <- sample_Y(I, K, Eta, Xi, na_indices, Yna)
    Xi          <- sample_Xi_dp(I, K, Eta, Xi, Y, com_siz, alpha_dp)
    #omega       <- sample_omega(K, alpha, Xi)
    sigsq       <- sample_sigsq(K, Q, a_sig, b_sig, U)
    omesq       <- sample_omesq(a_ome, b_ome, zeta)
    #no need to sample alpha right now
    #tmp         <- sample_alpha(b, n_tun_alpha, del2_alpha, n_alpha, n_burn, K, a_alpha, b_alpha, alpha, omega)
    #alpha       <- tmp$alpha; del2_alpha  <- tmp$del2_alpha; n_alpha     <- tmp$n_alpha; n_tun_alpha <- tmp$n_tun_alpha
    tmp         <- sample_U(b, n_tun_U, del2_U, n_U, n_burn, I, K, Q, sigsq, zeta, U, Xi, Yna)
    U           <- tmp$U
    del2_U      <- tmp$del2_U
    n_U         <- tmp$n_U
    n_tun_U     <- tmp$n_tun_U
    tmp         <- sample_zeta(b, n_tun_zeta, del2_zeta, n_zeta, n_burn, I, K, omesq, zeta, U, Xi, Yna)
    zeta        <- tmp$zeta
    del2_zeta   <- tmp$del2_zeta
    n_zeta      <- tmp$n_zeta
    n_tun_zeta  <- tmp$n_tun_zeta
    # ppp
    if ((b > n_burn) & (b%%n_skip == 0)) y_ppp <- y_ppp + Yna[na_indices]/n_sams
  }
  y_ppp
}

GOF <- function(Y, THETA, model, dataset, loc)  #no need to modify
{
  I <- get_I(Y)
  K <- (-1 + sqrt(1 + 8*ncol(THETA$Eta_chain)))/2
  B <- nrow(THETA$Eta_chain)
  gof <- WAIC(I, K, B, Y, THETA$Eta_chain, THETA$Xi_chain)
  save(gof, file = paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "gof.RData"))
}

incidence_matrix <- function(THETA) #no need to modify
{
  I <- ncol(THETA$Xi_chain)
  B <- nrow(THETA$Xi_chain)
  K <- (-1 + sqrt(1 + 8*ncol(THETA$Eta_chain)))/2
  A_vec <- incidence_matrix0(I, K, B, THETA$Xi_chain)
  A <- matrix(0, I, I)
  for (i in 1:(I-1)) for (ii in (i+1):I) A[i, ii] <- A_vec[get_k(i, ii, I)]
  A <- A + t(A)
  diag(A) <- 1
  ## plot
  pdf(paste0(loc, "outs/",dataset, "/", model, "_", dataset, "_", "incidence_matrix.pdf"), height = 4, width = 4, pointsize = 10)
  par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
  heat.plot0(mat = A, show.grid = F, tick = 0, labs = NA, main = "Incidence matrix")
  dev.off()
}

