mcmc_parallel <- function(Y, K_range, Q_range, n_sams, n_burn, n_skip, model, dataset, loc)
{
        ## info
        cat(" * MCMC", "\n", sep = "")
        ## indices
        KQ <- as.matrix(expand.grid(K_range, Q_range))
        ## set cores and cluster
        suppressMessages(suppressWarnings(library(doParallel)))
        cl <- makeCluster(min(detectCores(), length(KQ)))
        registerDoParallel(cl)
        tmp <- foreach (i = 1:nrow(KQ), .combine = rbind, .inorder = FALSE, .packages = c("Rcpp")) %dopar% {
                source   (paste0(loc, "code/", "rfunctions.R"))
                source   (paste0(loc, "code/", model, "_functions.R"))
                sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
                K <- KQ[i,1]
                Q <- KQ[i,2]
                MCMC(Y, K, Q, n_sams, n_burn, n_skip, model = paste0(model, "_K_", K, "_Q_", Q), dataset, loc)
        }
        stopCluster(cl)
}

gof_parallel <- function(Y, K_range, Q_range, model, dataset, loc)
{
        ## info
        cat(" * WAIC ", "\n", sep = "")
        ## indices
        KQ <- as.matrix(expand.grid(K_range, Q_range))
        ## set cores and cluster
        suppressMessages(suppressWarnings(library(doParallel)))
        cl <- makeCluster(min(detectCores(), nrow(KQ)))
        registerDoParallel(cl)  
        tmp <- foreach (i = 1:nrow(KQ), .combine = rbind, .inorder = FALSE, .packages = c("Rcpp")) %dopar% {
                source   (paste0(loc, "code/", "rfunctions.R"))
                source   (paste0(loc, "code/", model, "_functions.R"))
                sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
                K <- KQ[i,1]
                Q <- KQ[i,2]
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_", "samples.RData"))
                GOF(Y, THETA, paste0(model, "_K_", K, "_Q_", Q), dataset, loc)
        }
        stopCluster(cl)
}

cv_parallel <- function(Y, K_range, Q_range, n_burn, n_sams, n_skip, model, dataset, loc)
{
        ## info
        cat(" * CV", "\n", sep = "")
        ## indices
        L <- 5
        set.seed(1234)
        folds <- get_folds(M = nrow(Y), 1, L)
        KQL <- as.matrix(expand.grid(K_range, Q_range, 1:L))
        ## set cores and cluster
        suppressMessages(suppressWarnings(require(doParallel)))
        cl <- makeCluster(min(detectCores(), nrow(KQL)))
        registerDoParallel(cl)
        cv <- foreach (i = 1:nrow(KQL), .inorder = F, .packages = c("Rcpp", "mvtnorm", "gtools")) %dopar% {
                source   (paste0(loc, "code/", "rfunctions.R"))
                source   (paste0(loc, "code/", model, "_functions.R"))
                sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
                K <- KQL[i,1]
                Q <- KQL[i,2]
                l <- KQL[i,3]
                na_indices <- folds == l
                y_true <- as.factor(Y[na_indices])
                Yna <- Y
                Yna[na_indices] <- NA
                y_ppp <- YPPP(Yna, na_indices, K, Q, n_sams, n_burn, n_skip)
                cv <- list(l = l, K = K, y_true = y_true, y_ppp = y_ppp)
                save(cv, file = paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_l_", l, "_", dataset, "_", "cv.RData"))
        }
        stopCluster(cl)
}

test_parallel <- function(K_range, Q_range,model, dataset, loc)
{
        ## test statistics
        cat(" * TEST", "\n", sep = "")
        ## indices
        KQ <- as.matrix(expand.grid(K_range, Q_range))
        ## set cores and cluster
        suppressMessages(suppressWarnings(require(doParallel)))
        cl <- makeCluster(detectCores())
        registerDoParallel(cl)  
        test <- foreach (i = 1:nrow(KQ), .packages = c("Rcpp", "igraph"), .inorder = F) %dopar% {
                source   (paste0(loc, "code/rfunctions.R"))
                sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
                K <- KQ[i,1]
                Q <- KQ[i,2]
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_", "samples.RData"))
                I <- ncol(THETA$Xi_chain)
                B <- nrow(THETA$Xi_chain)
                stats <- matrix(NA, B, 6)
                for (b in 1:B) { 
                        Eta <- THETA$Eta_chain[b,]
                        Xi  <- THETA$Xi_chain[b,]
                        YY  <- simulate_data(I, K, Eta, Xi)
                        g   <- graph.adjacency(adjmatrix = YY, mode = "undirected", diag = F)
                        stats[b,] <- c(graph.density(g), transitivity(g), assortativity.degree(g), mean_distance(g), mean(degree(g)), sd(degree(g)))
                }
                save(K, Q, stats, file = paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_", "test.RData"))
        }
        stopCluster(cl)
}