mcmc_parallel <- function(Y, K_range, n_sams, n_burn, n_skip, model, dataset, loc)
{
        ## info
        cat(" * MCMC", "\n", sep = "")
        ## set cores and cluster
        suppressMessages(suppressWarnings(library(doParallel)))
        cl <- makeCluster(min(detectCores(), length(K_range)))
        registerDoParallel(cl)
        tmp <- foreach (i = 1:length(K_range), .combine = rbind, .inorder = FALSE, .packages = c("Rcpp")) %dopar% {
                source   (paste0(loc, "code/", "rfunctions.R"))
                source   (paste0(loc, "code/", model, "_functions.R"))
                sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
                K <- K_range[i]
                MCMC(Y, K, n_sams, n_burn, n_skip, model = paste0(model, "_K_", K), dataset, loc)
        }
        stopCluster(cl)
}

gof_parallel <- function(Y, K_range, model, dataset, loc)
{
        ## info
        cat(" * WAIC ", "\n", sep = "")
        ## set cores and cluster
        suppressMessages(suppressWarnings(library(doParallel)))
        cl <- makeCluster(min(detectCores(), length(K_range)))
        registerDoParallel(cl)  
        tmp <- foreach (i = 1:length(K_range), .combine = rbind, .inorder = FALSE, .packages = c("Rcpp")) %dopar% {
                source   (paste0(loc, "code/", "rfunctions.R"))
                source   (paste0(loc, "code/", model, "_functions.R"))
                sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
                K <- K_range[i]
                load(file = paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_", "samples.RData"))
                GOF(Y, THETA, paste0(model, "_K_", K), dataset, loc)
        }
        stopCluster(cl)
}

cv_parallel <- function(Y, K_range, n_burn, n_sams, n_skip, model, dataset, loc)
{
        ## info
        cat(" * CV ", "\n", sep = "")
        ## indices
        L <- 5
        set.seed(1234)
        folds <- get_folds(M = nrow(Y), 1, L)
        KL <- as.matrix(expand.grid(K_range, 1:L))
        ## set cores and cluster
        suppressMessages(suppressWarnings(require(doParallel)))
        cl <- makeCluster(min(detectCores(), nrow(KL)))
        registerDoParallel(cl)
        cv <- foreach (i = 1:nrow(KL), .inorder = F, .packages = c("Rcpp")) %dopar% {
                source   (paste0(loc, "code/", "rfunctions.R"))
                source   (paste0(loc, "code/", model, "_functions.R"))
                sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
                K <- KL[i,1]
                l <- KL[i,2]
                na_indices <- folds == l
                y_true <- as.factor(Y[na_indices])
                Yna <- Y
                Yna[na_indices] <- NA
                y_ppp <- YPPP(Yna, na_indices, K, n_sams, n_burn, n_skip)
                cv <- list(l = l, K = K, y_true = y_true, y_ppp = y_ppp)
                save(cv, file = paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_l_", l, "_", dataset, "_", "cv.RData"))
        }
        stopCluster(cl)
}

test_parallel <- function(K_range, model, dataset, loc)
{
        ## info
        cat(" * TEST ", "\n", sep = "")
        ## test statistics
        suppressMessages(suppressWarnings(require(doParallel)))
        cl <- makeCluster(detectCores())
        registerDoParallel(cl)  
        test <- foreach (i = 1:length(K_range), .packages = c("Rcpp", "igraph"), .inorder = F) %dopar% {
                source   (paste0(loc, "code/rfunctions.R"))
                sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
                K <- K_range[i]
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_", "samples.RData"))
                I <- ncol(THETA$Xi_chain)
                B <- nrow(THETA$Xi_chain)
                stats <- matrix(NA, B, 6)
                for (b in 1:B) { 
                        Lambda <- THETA$Lambda_chain[b,]
                        Xi     <- THETA$Xi_chain[b,]
                        YY     <- simulate_data(I, K, Lambda, Xi)
                        g      <- graph.adjacency(adjmatrix = YY, mode = "undirected", diag = F)
                        stats[b,] <- c(graph.density(g), transitivity(g), assortativity.degree(g), mean_distance(g), mean(degree(g)), sd(degree(g)))
                }
                save(K, stats, file = paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_", "test.RData"))
        }
        stopCluster(cl)
}