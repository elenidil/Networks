generate_data <- function(I, K, omega, mu_mu, sig2_mu, a_sig, b_sig, model, dataset, loc, plot_data = T)
{
        suppressMessages(suppressWarnings(library(igraph)))
        ## parameters
        Xi     <- sort(sample(x = 1:K, size = I, replace = T, prob = omega))
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
        omega  <<- omega
        mu     <<- mu
        sigsq  <<- sigsq
        Lambda <<- Lambda
        Y      <<- Y
        Ycube  <<- Ycube
        ll     <<- loglik(I, K, Lambda, Xi, Y)
        ## return
        list(Y = Y, Ycube = Ycube, omega = omega, Xi = Xi, mu = mu, sigsq = sigsq, Lambda = Lambda, ll = ll)
}