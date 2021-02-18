generate_data <- function(I, K, Q, sigsq, zeta, Xi, plot_data = T)
{
        ## parameters
        suppressMessages(suppressWarnings(library(mvtnorm)))
        suppressMessages(suppressWarnings(library(igraph)))
        U   <- rmvnorm(n = K, mean = rep(0, Q), sigma = diag(sigsq, Q))
        Eta <- matrix(NA, c(K*(K+1)/2), 1)
        for(k in 1:K) {
                for(kk in k:K) {
                        Eta[get_k_diag(k, kk, K)] <- zeta - sqrt(sum((U[k,] - U[kk,])^2))
                }
        }
        ## Likelihood
        iPcube <- Ycube <- array(0, c(I, I))
        for(i in 1:(I-1)) {
                for(ii in (i+1):I) {
                        iPcube[i, ii] <- iPcube[ii, i] <- expit(Eta[get_k_diag(min(Xi[i], Xi[ii]), max(Xi[i], Xi[ii]), K)])
                        Ycube [i, ii] <- Ycube [ii, i] <- rbinom(n = 1, size = 1, prob = iPcube[i, ii])
                }
        }
        iP <- Y <- matrix(NA, I*(I-1)/2, 1)
        for (i in 1:(I-1)) {
                for (ii in (i+1):I) {
                        iP[get_k(i, ii, I)] <- iPcube[i, ii]
                        Y [get_k(i, ii, I)] <- Ycube [i, ii]
                }
        }
        ## plot
        if (plot_data) {
                pdf(paste0(loc, "outs/", model, "_", dataset,"_adjacency.pdf"), height = 1*4, width = 1*4, pointsize = 10)
                par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
                adjacency.plot(mat = Ycube,  labs = NA, show.grid = F, tick = 0, main = "Adjacency Matrix")
                dev.off()
                
                pdf(paste0(loc, "outs/", model, "_", dataset,"_probabilities.pdf"), height = 1*4, width = 1*4, pointsize = 10)
                par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
                heat.plot0    (mat = iPcube, labs = NA, show.grid = F, tick = 0, main = "Interaction Probabilities")
                dev.off()
                
                pdf(paste0(loc, "outs/", model, "_", dataset,"_graph.pdf"), height = 1*4, width = 1*4, pointsize = 10)
                par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
                g <- graph.adjacency(adjmatrix = Ycube, mode = "undirected", diag = F)
                plot(g, vertex.color = "lightblue", edge.color = "gray50", vertex.label = NA, vertex.size = 3)
                title("Graph")
                box()
                dev.off()        
        }
        ## load
        Y      <<- Y
        Ycube  <<- Ycube
        iPcube <<- iPcube
        U      <<- U
        Xi     <-  Xi-1
        Xi     <<- Xi
        Eta    <<- Eta
        zeta   <<- zeta
        sigsq  <<- sigsq
        ll     <<- loglik(I, K, Eta, Xi, Y)
        ## return
        list(Y = Y, Ycube = Ycube, Xi = Xi, U = U, Eta = Eta, zeta = zeta, sigsq = sigsq, ll = ll)
}