################################### PLOTS ######################################
adjacency.plot <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis, ...)
{ 
        JJ <- dim(mat)[1]
        colorscale <- c("white", rev(heat.colors(100)))
        if(missing(labs))     labs <- 1:JJ
        if(missing(col.axis)) col.axis <- rep("black", JJ)
        if(missing(cex.axis)) cex.axis <- 0.5
        if(missing(tick))     tick <- TRUE
        ## adjacency matrix
        image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
              col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
        for(j in 1:JJ){
                axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
                axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
        }
        box()
        if(show.grid) grid(nx=JJ, ny=JJ)
}

heat.plot0 <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis,...)
{ 
        JJ <- dim(mat)[1]
        colorscale <- c("white", rev(heat.colors(100)))
        if(missing(labs))     labs <- 1:JJ
        if(missing(col.axis)) col.axis <- rep("black", JJ)
        if(missing(cex.axis)) cex.axis <- 0.5
        if(missing(tick))     tick <- TRUE
        ## adjacency matrix
        image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
              col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
        for(j in 1:JJ){
                axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
                axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
        }
        box()
        if(show.grid) grid(nx = JJ, ny = JJ)
}

plot_chain0 <- function(chain, plot_ch = T, subset_indices = F, ...)
{
        chain <- as.matrix(chain)
        npars <- ncol(chain)
        npars <- ncol(chain)
        if (subset_indices) {
                n <- 6
                if (n > npars) 
                        indices <- 1:npars
                else 
                        indices <- sort(sample(x = 1:npars, size = n, replace = F))
                param <- as.matrix(c(param)[indices])
                chain <- as.matrix(chain[,indices])
                npars <- ncol(chain)
        }
        if (plot_ch) {
                nr <- ceiling(npars/2)
                nc <- 2 - as.numeric(npars==1)
                windows(width = 10*nc, height = 5*nr)
                par(mfrow = c(nr, nc))
                for (i in 1:npars) {
                        plot(chain[,i], type = "l", col = "gray",
                             cex.axis = 0.7, cex.lab = 0.7,
                             ylim = range(chain[,i]),
                             xlab = "Iteration", ylab = "Parameter", ...)
                        abline(h = mean(chain[,i]), col = 1, lwd = 2, lty = 1)
                        abline(h = quantile(chain[,i], c(0.025, 0.975)), col = 1, lwd = 1, lty = 2)
                }
        }
}

plot_chain <- function(chain, param, plot_ch = T, subset_indices = F, ...)
{
        param <- as.matrix(c(param))
        chain <- as.matrix(chain)
        npars <- ncol(chain)
        npars <- ncol(chain)
        if (subset_indices) {
                n <- 6
                if (n > npars) 
                        indices <- 1:npars
                else 
                        indices <- sort(sample(x = 1:npars, size = n, replace = F))
                param <- as.matrix(c(param)[indices])
                chain <- as.matrix(chain[,indices])
                npars <- ncol(chain)
        }
        if (plot_ch) {
                nr <- ceiling(npars/2)
                nc <- 2 - as.numeric(npars==1)
                windows(width = 10*nc, height = 5*nr)
                par(mfrow = c(nr, nc))
                for (i in 1:npars) {
                        plot(chain[,i], type = "l", col = "gray",
                             cex.axis = 0.7, cex.lab = 0.7,
                             ylim = range(chain[,i], param[i]),
                             xlab = "Iteration", ylab = "Parameter", ...)
                        abline(h = mean(chain[,i]), col = 1, lwd = 2, lty = 1)
                        abline(h = quantile(chain[,i], c(0.025, 0.975)), col = 1, lwd = 1, lty = 2)
                        abline(h = param[i], col = 2, lwd = 2, lty = 3)
                }
        }
}

plot_ci <- function(chain, param, subset_indices = F, n = 25)
{
        param <- as.matrix(c(param))
        chain <- as.matrix(chain)
        if (subset_indices) {
                if (missing(n)) n <- 25
                indices <- sort(sample(1:ncol(chain), size = n))
                param   <- as.matrix(param[indices,])
                chain   <- as.matrix(chain[,indices]) 
        }
        npars <- ncol(chain)

        PM  <- as.matrix(colMeans(chain))
        CI1 <- t(apply(X = chain, MARGIN = 2, FUN = function(x) quantile(x, probs = c(0.025, 0.975))))
        CI2 <- t(apply(X = chain, MARGIN = 2, FUN = function(x) quantile(x, probs = c(0.005, 0.995))))
        
        windows()
        xrange <- c(1, npars)
        yrange <- range(CI2)
        plot(NA, NA, type = "n", xlim = xrange, ylim = yrange, xaxt = "n", xlab = "Parameter", ylab = "Parameter value")
        segments(x0 = 1:npars, y0 = CI1[,1], x1 = 1:npars, y1 = CI1[,2], lwd = 2)
        segments(x0 = 1:npars, y0 = CI2[,1], x1 = 1:npars, y1 = CI2[,2], lwd = 1)
        lines(x = 1:npars, y = PM,    type = "p", pch = 16, cex = 0.7, col = 1)
        lines(x = 1:npars, y = param, type = "p", pch = 18, cex = 0.7, col = 2)
}

get_measures_0 <- function(Ycube, model, dataset, loc)
{
        suppressMessages(suppressWarnings(require(PRROC)))
        suppressMessages(suppressWarnings(library(igraph)))
        ## sizes
        L <- 5
        ## observed statistics
        g <- graph.adjacency(adjmatrix = Ycube, mode = "undirected", diag = F)
        stats_obs <- c(graph.density(g), transitivity(g), assortativity.degree(g), mean_distance(g), mean(degree(g)), sd(degree(g)))
        n_stats <- length(stats_obs)
        ## table
        nam <- c("K", "Q", "AUC", "WAIC1", "WAIC2", "Density", "Transitivity", "Assortativity", "Distance", "Mean Degree", "SD Degree")
        tab <- as.data.frame(matrix(NA, 1, length(nam)))
        colnames(tab) <- nam
        ## auc
        tmp <- matrix(NA, L, 2)
        for (l in 1:L) {
                load(paste0(loc, "outs/", dataset, "/", model, "_l_", l, "_", dataset, "_", "cv.RData"))
                tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
        }
        tab[3] <- round(mean(tmp[,1]), 6)
        ## gof
        load(paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "gof.RData"))
        tab[4] <- round(gof$waic1, 3)
        tab[5] <- round(gof$waic2, 3)
        ## test
        load(paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "test.RData"))
        for (l in 1:n_stats) tab[5+l] <- round(mean(stats[,l] > stats_obs[l], na.rm = T), 6)
        ## save
        save(tab, file = paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "measures.RData"))
}

get_measures_1 <- function(Ycube, K_range, model, dataset, loc)
{
        suppressMessages(suppressWarnings(require(PRROC)))
        suppressMessages(suppressWarnings(library(igraph)))
        L <- 5
        n_K <- length(K_range)
        ## observed statistics
        g <- graph.adjacency(adjmatrix = Ycube, mode = "undirected", diag = F)
        stats_obs <- c(graph.density(g), transitivity(g), assortativity.degree(g), mean_distance(g), mean(degree(g)), sd(degree(g)))
        n_stats <- length(stats_obs)
        ## table
        nam <- c("K", "Q", "AUC", "WAIC1", "WAIC2", "Density", "Transitivity", "Assortativity", "Distance", "Mean Degree", "SD Degree")
        tab <- as.data.frame(matrix(NA, n_K, length(nam)))
        colnames(tab) <- nam
        for (i in 1:n_K) {
                ## dim
                tab[i, 1] <- K <- K_range[i]
                ## auc
                tmp <- NULL
                for (l in 1:L) {
                        load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_l_", l, "_", dataset, "_", "cv.RData"))
                        tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
                }
                tab[i, 3] <- round(mean(tmp), 6)
                ## gof
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_", "gof.RData"))
                tab[i, 4] <- round(gof$waic1, 3)
                tab[i, 5] <- round(gof$waic2, 3)
                ## test
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_", "test.RData"))
                for (l in 1:n_stats) tab[i, 5+l] <- round(mean(stats[,l] > stats_obs[l], na.rm = T), 6)
        }
        save(tab, file = paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "measures.RData"))
}

get_measures_2 <- function(Ycube, K_range, Q_range, model, dataset, loc)
{
        suppressMessages(suppressWarnings(require(PRROC)))
        suppressMessages(suppressWarnings(library(igraph)))
        ## sizes
        L <- 5
        n_K <- length(K_range)
        n_Q <- length(Q_range)
        ## observed statistics
        g <- graph.adjacency(adjmatrix = Ycube, mode = "undirected", diag = F)
        stats_obs <- c(graph.density(g), transitivity(g), assortativity.degree(g), mean_distance(g), mean(degree(g)), sd(degree(g)))
        n_stats <- length(stats_obs)
        ## table
        QK <- as.matrix(expand.grid(Q_range, K_range))
        nam <- c("K", "Q", "AUC", "WAIC1", "WAIC2", "Density", "Transitivity", "Assortativity", "Distance", "Mean Degree", "SD Degree")
        tab <- as.data.frame(matrix(NA, nrow(QK), length(nam)))
        colnames(tab) <- nam
        for (i in 1:nrow(QK)) {
                tab[i, 1] <- K <- QK[i,2]
                tab[i, 2] <- Q <- QK[i,1]
                ## auc
                tmp <- NULL
                for (l in 1:L) {
                        load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_l_", l, "_", dataset, "_", "cv.RData"))
                        tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
                }
                tab[i, 3] <- round(mean(tmp), 6)
                ## gof
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_", "gof.RData"))
                tab[i, 4] <- round(gof$waic1, 2)
                tab[i, 5] <- round(gof$waic2, 2)
                ## test
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_", "test.RData"))
                for (l in 1:n_stats) tab[i, 5+l] <- round(mean(stats[,l] > stats_obs[l], na.rm = T), 6)
        }
        save(tab, file = paste0(loc, "outs/", dataset, "/", model, "_", dataset, "_", "measures.RData"))
}

################################## GENERIC #######################################

dec <- function(x, k) format(format(round(x, k), nsmall = k))

expit <- function(x) 
{
        1/(1 + exp(-x))
}

procus <- function(Z, Z0, K)
{
        ## Z closest to Z0
        for(i in 1:K) Z[,i] <- Z[,i]-mean(Z[,i]) + mean(Z0[,i])  # translation
        A     <- t(Z)%*%(Z0%*%t(Z0))%*%Z
        eA    <- eigen(A, symmetric=T)
        Ahalf <- eA$vec[,1:K]%*%diag(sqrt(eA$val[1:K]))%*%t(eA$vec[,1:K])
        t(t(Z0)%*%Z%*%solve(Ahalf)%*%t(Z)) 
}

################################# NETWORKS ##########################################

get_I <- function(Y) (1 + sqrt(1 + 8*nrow(Y)))/2

get_k <- function(i ,ii, I)
{
        I*(I-1)/2 - (I-i+1)*(I-i)/2 + ii - i 
}        

get_k_diag <- function(k, kk, K)
{
        K*(k-1) + kk - (k-1)*k/2  
}

get_folds <- function(M, J, L)
{
        folds <- matrix(NA, M, J)
        for (j in 1:J) {
                for (l in 1:L) {
                        indices <- which(is.na(folds[,j]))
                        folds[sample(indices, floor(M/L), F), j] <- l
                }
                indices <- which(is.na(folds[,j]))
                if (sum(indices)) folds[indices, j] <- sample(1:L, 1, F)
        }
        folds
}

latent_colors <- function(A)
{
        r <- atan2(A[,2], A[,1])
        r <- r + abs(min(r))
        r <- r/max(r)
        g <- 1 - r
        b <- A[,2]^2 + A[,1]^2
        b <- b/max(b)
        list(r = r, g = g, b = b)
}