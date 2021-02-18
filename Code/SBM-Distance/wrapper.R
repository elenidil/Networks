wrapper_0 <- function(n_sams, n_burn, n_skip, model, dataset, loc)
{
        ## info
        cat("Model   : ", model,   "\n", "Dataset : ", dataset, "\n", sep = "")
        
        ## LIBS
        suppressMessages(suppressWarnings(library(Rcpp)))
        
        ## PATHS AND SOURCE 
        source   (paste0(loc, "code/", "path.R"))
        source   (paste0(loc, "code/", "rfunctions.R"))
        source   (paste0(loc, "code/", model, "_parallel.R"))
        source   (paste0(loc, "code/", model, "_functions.R"))
        sourceCpp(paste0(loc, "code/", model, "_functions.cpp"))
        
        ## DATA
        load(paste0(loc, "data/", dataset, "_data.RData"))
        
        ## MCMC
        MCMC(Y, n_sams, n_burn, n_skip, model, dataset, loc)
        
        ## GOODNESS-OF-FIT
        GOF(Y, model, dataset, loc)
        
        ## CROSS-VALIDATION
        cv_parallel(Y, n_burn, n_sams, n_skip, model, dataset, loc)
        
        ## POSTERIOR PREDICTIVE P-VALUES
        test_parallel(I, model, dataset, loc)
        
        ## RESUTLS
        get_measures_0(Ycube, model, dataset, loc)
        
        TRUE
}

wrapper_1 <- function(K_range, n_sams, n_burn, n_skip, model, dataset, loc)
{
        ## info
        cat("Model   : ", model,   "\n", "Dataset : ", dataset, "\n", sep = "")
        
        ## PATHS AND SOURCE 
        source(paste0(loc, "code/", model, "_parallel.R"))
        source(paste0(loc, "code/", "path.R"))
        source(paste0(loc, "code/", "rfunctions.R"))

        ## DATA
        load(paste0(loc, "data/", dataset, "_data.RData"))
        
        ## MCMC
        mcmc_parallel(Y, K_range, n_sams, n_burn, n_skip, model, dataset, loc)
        
        ## GOODNESS-OF-FIT
        gof_parallel(Y, K_range, model, dataset, loc)
        
        ## CROSS-VALIDATION
        cv_parallel(Y, K_range, n_burn, n_sams, n_skip, model, dataset, loc)
        
        ## POSTERIOR PREDICTIVE P-VALUES
        test_parallel(K_range, model, dataset, loc)
        
        ## RESUTLS
        get_measures_1(Ycube, K_range, model, dataset, loc)
        
        TRUE
}

wrapper_2 <- function(K_range, Q_range, n_sams, n_burn, n_skip, model, dataset, loc)
{
        ## info
        cat("Model   : ", model,   "\n", "Dataset : ", dataset, "\n", sep = "")
        
        ## PATHS AND SOURCE 
        source(paste0(loc, "code/", model, "_parallel.R"))
        source(paste0(loc, "code/", "path.R"))
        source(paste0(loc, "code/", "rfunctions.R"))

        ## DATA
        load(paste0(loc, "data/", dataset, "_data.RData"))
        
        ## MCMC
        mcmc_parallel(Y, K_range, Q_range, n_sams, n_burn, n_skip, model, dataset, loc)
        
        ## GOODNESS-OF-FIT
        gof_parallel(Y, K_range, Q_range, model, dataset, loc)
        
        ## CROSS-VALIDATION
        cv_parallel(Y, K_range, Q_range, n_burn, n_sams, n_skip, model, dataset, loc)
        
        ## POSTERIOR PREDICTIVE P-VALUES
        test_parallel(K_range, Q_range, model, dataset, loc)
        
        ## RESUTLS
        get_measures_2(Ycube, K_range, Q_range, model, dataset, loc)
        
        TRUE
}

get_results <- function(dataset, loc)
{
        out <- NULL
        models <- c("erdos", "dist", "eigen", "class", "clist", "cleigen", "class2")
        for (i in 1:length(models)) {
                tfile <- paste0(loc, "outs/", dataset, "/", models[i], "_", dataset, "_", "measures.RData") 
                if (file.exists(tfile)) load(tfile)
                out <- rbind(out, data.frame(model = models[i], tab))
                rm(tab)
        }
        out <- as.data.frame(out)
        colnames(out) <- c("model", "K", "Q", "AUC", "WAIC1", "WAIC2", "Density", "Transitivity", "Assortativity", "Distance", "Mean Degree", "SD Degree")
        tab <- out
        save(tab, file = paste0(loc, "outs/", dataset, "/", dataset, "_", "measures.RData"))
        write.table(x = tab, file = paste0(loc, "outs/", dataset, "/", dataset, "_", "measures.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
        tab
}

plot_measures <- function(K_range, Q_range, dataset, loc)
{
        suppressMessages(suppressWarnings(require(PRROC)))
        L   <- 5
        n_K <- length(K_range)
        n_Q <- length(Q_range)
        K_models <- c("eigen", "dist", "class")
        Q_models <- c("clist", "cleigen", "class2")
        ## erdos
        tmp <- NULL
        for (l in 1:L) {
                load(paste0(loc, "outs/", dataset, "/", "erdos", "_l_", l, "_", dataset, "_", "cv.RData"))
                tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
        }
        y_0 <- mean(tmp)
        ## latent models
        y_1 <- NULL
        for (h in 1:length(K_models)) {
                auc <- NULL
                for (i in 1:n_K) {
                        model <- K_models[h]
                        K <- K_range[i]
                        tmp <- NULL
                        for (l in 1:L) {
                                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_l_", l, "_", dataset, "_", "cv.RData"))
                                tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
                        }
                        auc[i] <- mean(tmp)
                }
                y_1[h] <- max(auc)
        }
        ## multilevel models
        y_2 <- NULL
        for (h in 1:length(Q_models)) {
                auc <- NULL
                for (i in 1:n_K) {
                        for (j in 1:n_Q) {
                                model <- Q_models[h]
                                K <- K_range[i]
                                Q <- Q_range[j]
                                tmp <- NULL
                                for (l in 1:L) {
                                        load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_l_", l, "_", dataset, "_", "cv.RData"))
                                        tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
                                }
                                auc <- mean(tmp)
                        }
                }
                y_2[h] <- max(auc)
        }
        ## plot
        y <- c(y_1, y_2, y_0)
        models <- c(K_models, Q_models, "erdos")
        n_y <- length(y)
        mycolors <- c("dodgerblue2", "red3", "green3", "slateblue", "darkorange", "skyblue1", "violetred4")
        mypch <- c(0,1,2,5,6,3,4)
        pdf(paste0(loc, "outs/", dataset, "/", dataset, "_", "measures.pdf"), height = 6, width = 6, pointsize = 12)
        par(mfrow = c(1, 1), mar = c(3, 4, 3, 2) - 0, mgp = c(2, 1, 0), oma = 0.1*c(1, 1, 1, 1))
        plot(1:n_y, y, type = "p", lty = 1, pch = mypch[1:n_y], col = mycolors[1:n_y], xlim = c(0.5,n_y+0.5),
             ylim = c(0.45, 0.975), xaxt = "n", xlab = "Model", ylab = "AUC", cex = 1.1, cex.axis = 0.7)
        axis(side = 1, at = 1:n_y, labels = models, cex.axis = 0.7)
        text(x = 1:n_y, y = y, labels = round(y, 3), pos = 3, col = mycolors[1:n_y], cex = 0.8)
        legend("bottomleft", legend = models, lty = 1, col = mycolors[1:n_y], pch = mypch[1:n_y], ncol = 1, bty = "n")
        title(main = dataset, outer = T, line = -1)
        dev.off()
}

plot_measures2 <- function(K_range, Q_range, dataset, loc)
{
        suppressMessages(suppressWarnings(require(PRROC)))
        suppressMessages(suppressWarnings(library(ggplot2)))
        L   <- 5
        n_K <- length(K_range)
        n_Q <- length(Q_range)
        K_models <- c("eigen", "dist", "class")
        Q_models <- c("clist", "cleigen", "class2")
        ## erdos
        tmp <- NULL
        for (l in 1:L) {
                load(paste0(loc, "outs/", dataset, "/", "erdos", "_l_", l, "_", dataset, "_", "cv.RData"))
                tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
        }
        y_0 <- as.matrix(tmp)
        ## latent models
        y_1 <- NULL
        for (h in 1:length(K_models)) {
                auc <- NULL
                for (i in 1:n_K) {
                        model <- K_models[h]
                        K <- K_range[i]
                        tmp <- NULL
                        for (l in 1:L) {
                                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_l_", l, "_", dataset, "_", "cv.RData"))
                                tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
                        }
                        auc <- c(auc, tmp)
                }
                y_1 <- cbind(y_1, auc)
        }
        ## multilevel models
        y_2 <- NULL
        for (h in 1:length(Q_models)) {
                auc <- NULL
                for (i in 1:n_K) {
                        for (j in 1:n_Q) {
                                model <- Q_models[h]
                                K <- K_range[i]
                                Q <- Q_range[j]
                                tmp <- NULL
                                for (l in 1:L) {
                                        load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_l_", l, "_", dataset, "_", "cv.RData"))
                                        tmp[l] <- roc.curve(scores.class0 = cv$y_ppp, weights.class0 = as.numeric(cv$y_true)-1)$auc
                                }
                                auc <- c(auc, tmp)
                        }
                }
                y_2 <- cbind(y_2, auc)
        }
        ## plot
        models <- c("erdos", "eigen", "dist", "class", "c-dist", "c-eigen", "c-class")
        auc <- data.frame(x = c(rep(1, nrow(y_0)), rep(2:4, rep(nrow(y_1), 3)), rep(5:7, rep(nrow(y_2), 3))), auc = c(c(y_0), c(y_1), c(y_2)))
        jitter <- position_jitter(width = 0.1)
        ggplot(data = auc, mapping = aes(x, auc)) + 
                geom_point(size = 1, alpha = 0.5, position = jitter) + 
                labs(color = "Model", x = "Model", y = "AUC") + 
                scale_x_discrete(limits = models) +
                ylim(0, 1)
        ggsave(filename = paste0(loc, "outs/", dataset, "/", dataset, "_", "measures_2.pdf"), height = 6, width = 6)
}

plot_stats2 <- function(K_range, Q_range, dataset, loc)
{
        suppressMessages(suppressWarnings(library(ggplot2)))
        suppressMessages(suppressWarnings(library(gridExtra)))
        suppressMessages(suppressWarnings(library(igraph)))
        n_K <- length(K_range)
        n_Q <- length(Q_range)
        K_models <- c("eigen", "dist", "class")
        Q_models <- c("clist", "cleigen", "class2")
        K_models_disp <- c("eigen", "dist", "class")
        Q_models_disp <- c("c-dist", "c-eigen", "c-class")
        features <- c("density", "transitivity", "assortativity", "mean distance", "mean degree", "sd degree")
        models   <- c("erdos", K_models_disp, Q_models_disp)
        ## data
        load(paste0(loc, "data/", dataset, "_data.RData"))
        ## observed statistics
        g <- graph.adjacency(adjmatrix = Ycube, mode = "undirected", diag = F)
        stats_obs <- c(graph.density(g), transitivity(g), assortativity.degree(g), mean_distance(g), mean(degree(g)), sd(degree(g)))
        n_stats <- length(stats_obs)
        ## stats
        ic <- NULL
        for (h in 1:n_stats) {
                feat <- features[h]
                # erdos
                model <- "erdos"
                load(paste0(loc, "outs/", dataset, "/", "erdos_", dataset, "_", "test.RData"))
                ic <- rbind(ic, c(feat, model, "95",  quantile(stats[,h], c(0.025, 0.975)), mean(stats[,h]), stats_obs[h]))
                ic <- rbind(ic, c(feat, model, "99",  quantile(stats[,h], c(0.005, 0.995)), mean(stats[,h]), stats_obs[h]))
                rm(model, stats)
                # K models
                for (k in 1:length(K_models)) {
                        model <- K_models[k]
                        waic  <- NULL
                        for (i in 1:n_K) {
                                K <- K_range[i]
                                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_gof.RData"))
                                waic <- c(waic, gof$waic1) 
                                rm(K, gof)
                        }
                        K <- K_range[(!is.nan(waic)) & (waic == min(waic, na.rm = T))]
                        load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_test.RData"))
                        ic <- rbind(ic, c(feat, K_models_disp[k], "95",  quantile(stats[,h], c(0.025, 0.975)), mean(stats[,h]), stats_obs[h]))
                        ic <- rbind(ic, c(feat, K_models_disp[k], "99",  quantile(stats[,h], c(0.005, 0.995)), mean(stats[,h]), stats_obs[h]))
                        rm(model, K, stats)
                }
                # Q models
                for (q in 1:length(Q_models)) {
                        model <- Q_models[q]
                        tmp  <- matrix(NA, n_K*n_Q, 3)
                        for (i in 1:n_K) {
                                for (j in 1:n_Q) {
                                        K <- K_range[i]
                                        Q <- K_range[j]
                                        load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_gof.RData"))
                                        tmp[(i-1)*n_K+j,] <- c(K, Q, gof$waic1)
                                        rm(K, Q, gof)
                                }
                        }
                        waic <- tmp[,3]
                        K <- tmp[which((!is.nan(waic)) & (waic == min(waic, na.rm = T))), 1]
                        Q <- tmp[which((!is.nan(waic)) & (waic == min(waic, na.rm = T))), 2]
                        load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_test.RData"))
                        ic <- rbind(ic, c(feat, Q_models_disp[q], "95",  quantile(stats[,h], c(0.025, 0.975)), mean(stats[,h]), stats_obs[h]))
                        ic <- rbind(ic, c(feat, Q_models_disp[q], "99",  quantile(stats[,h], c(0.005, 0.995)), mean(stats[,h]), stats_obs[h]))
                        rm(model, K, Q, stats)
                }
                rm(feat)
        }
        ic <- as.data.frame(ic)
        colnames(ic) <- c("stat", "model", "conf", "lower", "upper", "mean", "obs")
        ic$lower <- as.numeric(as.character(ic$lower))
        ic$upper <- as.numeric(as.character(ic$upper))
        ic$mean  <- as.numeric(as.character(ic$mean))
        ic$obs   <- as.numeric(as.character(ic$obs))
        ## plots
        # for (h in 1:n_stats) {
        #         ic_h_95 <- ic[(ic$stat == features[h]) & (ic$conf == "95"), c("model", "lower", "upper", "mean", "obs")]
        #         ic_h_99 <- ic[(ic$stat == features[h]) & (ic$conf == "99"), c("model", "lower", "upper", "mean", "obs")]
        #         ggplot() + 
        #                 geom_pointrange(data = ic_h_95, mapping = aes(x = factor(model, levels = model), y = mean, ymin = upper, ymax = lower, ), shape = 15, fatten = 1, size = 0.5) +
        #                 geom_pointrange(data = ic_h_99, mapping = aes(x = model, y = mean, ymin = lower, ymax = upper), shape = 15, fatten = 1, size = 0.25) + 
        #                 geom_point(data = ic_h_95, mapping = aes(x = model, y = obs), shape = 16, size = 1.75, col = "royalblue") +
        #                 labs(x = "Model", y = "Credible Interval")
        #         ggsave(filename = paste0(loc, "outs/", dataset, "/", dataset, "_", "stat_", features[h], "_2.pdf"), height = 6, width = 6)
        # }
        for (h in 1:n_stats) {
                ic_h_95 <- ic[(ic$stat == features[h]) & (ic$conf == "95"), c("model", "lower", "upper", "mean", "obs")]
                ic_h_99 <- ic[(ic$stat == features[h]) & (ic$conf == "99"), c("model", "lower", "upper", "mean", "obs")]
                assign(paste0("p", h), 
                       ggplot() + 
                               geom_pointrange(data = ic_h_95, mapping = aes(x = factor(model, levels = model), y = mean, ymin = upper, ymax = lower, ), shape = 15, fatten = 1, size = 0.5) +
                               geom_pointrange(data = ic_h_99, mapping = aes(x = model, y = mean, ymin = lower, ymax = upper), shape = 15, fatten = 1, size = 0.25) + 
                               geom_point(data = ic_h_95, mapping = aes(x = model, y = obs), shape = 16, size = 1.75, col = "royalblue") +
                               labs(x = "Model", y = "Credible Interval", title = features[h])) + theme(text = element_text(size = rel(1)))
        }          
        ggsave(filename = paste0(loc, "outs/", dataset, "/", dataset, "_", "stat_2.pdf"), plot = grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2), height = 8, width = 12)
}

plot_rand2 <- function(K_range, Q_range, dataset, loc, B)
{
        # https://davetang.org/muse/2017/09/21/adjusted-rand-index/
        # load package
        suppressMessages(suppressWarnings(library(clues)))
        suppressMessages(suppressWarnings(library(mltools)))
        suppressMessages(suppressWarnings(library(ggplot2)))
        suppressMessages(suppressWarnings(library(gridExtra)))
        suppressMessages(suppressWarnings(library(Rcpp)))
        sourceCpp(paste0(loc, "code/get_incidence_matrix.cpp"))
        n_K <- length(K_range)
        n_Q <- length(Q_range)
        K_models <- c("eigen", "dist", "class")
        Q_models <- c("clist", "cleigen", "class2")
        K_models_disp <- c("eigen", "dist", "class")
        Q_models_disp <- c("c-dist", "c-eigen", "c-class")
        features <- c("Rand", "HA", "MA", "FM", "Jaccard", "MCC")
        models   <- c("erdos", K_models_disp, Q_models_disp)
        # ground truth
        load(paste0(loc, "data/", dataset, "_labels.RData"))
        gti <- c(get_incidence_matrix(gt))
        ic <- NULL
        # class model
        model <- "class"
        waic  <- NULL
        for (i in 1:n_K) {
                K <- K_range[i]
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_gof.RData"))
                waic <- c(waic, gof$waic1) 
                rm(K, gof)
        }
        K <- K_range[(!is.nan(waic)) & (waic == min(waic, na.rm = T))]
        load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_", dataset, "_samples.RData"))
        tmp <- NULL
        for (b in 1:B) tmp <- rbind(tmp, c(adjustedRand(cl1 = gt, cl2 = THETA$Xi_chain[b,]), mcc(c(get_incidence_matrix(THETA$Xi_chain[b,])), gti)))
        ic <- rbind(ic, cbind(features, model, "95", t(apply(tmp, 2, function(x) c(quantile(x, probs = c(0.025, 0.975)), mean(x))))))
        ic <- rbind(ic, cbind(features, model, "99", t(apply(tmp, 2, function(x) c(quantile(x, probs = c(0.005, 0.995)), mean(x))))))
        rm(model, K, THETA)
        # Q models
        for (q in 1:length(Q_models)) {
                model <- Q_models[q]
                tmp  <- matrix(NA, n_K*n_Q, 3)
                for (i in 1:n_K) {
                        for (j in 1:n_Q) {
                                K <- K_range[i]
                                Q <- K_range[j]
                                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_gof.RData"))
                                tmp[(i-1)*n_K+j,] <- c(K, Q, gof$waic1)
                                rm(K, Q, gof)
                        }
                }
                waic <- tmp[,3]
                K <- tmp[which((!is.nan(waic)) & (waic == min(waic, na.rm = T))), 1]
                Q <- tmp[which((!is.nan(waic)) & (waic == min(waic, na.rm = T))), 2]
                load(paste0(loc, "outs/", dataset, "/", model, "_K_", K, "_Q_", Q, "_", dataset, "_samples.RData"))
                tmp <- NULL
                for (b in 1:B) tmp <- rbind(tmp, c(adjustedRand(cl1 = gt, cl2 = THETA$Xi_chain[b,]), mcc(c(get_incidence_matrix(THETA$Xi_chain[b,])), gti)))
                ic <- rbind(ic, cbind(features, Q_models_disp[q], "95", t(apply(tmp, 2, function(x) c(quantile(x, probs = c(0.025, 0.975)), mean(x))))))
                ic <- rbind(ic, cbind(features, Q_models_disp[q], "99", t(apply(tmp, 2, function(x) c(quantile(x, probs = c(0.005, 0.995)), mean(x))))))
                rm(K, Q, THETA)
        }
        ic <- as.data.frame(ic)
        colnames(ic) <- c("stat", "model", "conf", "lower", "upper", "mean")
        rownames(ic) <- NULL
        ic$lower <- as.numeric(as.character(ic$lower))
        ic$upper <- as.numeric(as.character(ic$upper))
        ic$mean  <- as.numeric(as.character(ic$mean))
        for (h in 1:length(features)) {
                ic_h_95 <- ic[(ic$stat == features[h]) & (ic$conf == "95"), c("model", "lower", "upper", "mean")]
                ic_h_99 <- ic[(ic$stat == features[h]) & (ic$conf == "99"), c("model", "lower", "upper", "mean")]
                assign(paste0("p", h), 
                       ggplot() + 
                               geom_pointrange(data = ic_h_95, mapping = aes(x = factor(model, levels = model), y = mean, ymin = upper, ymax = lower, ), shape = 15, fatten = 1, size = 0.5) +
                               geom_pointrange(data = ic_h_99, mapping = aes(x = model, y = mean, ymin = lower, ymax = upper), shape = 15, fatten = 1, size = 0.25) + 
                               labs(x = "Model", y = "Credible Interval", title = features[h])) + theme(text = element_text(size = rel(1)))
        }          
        ggsave(filename = paste0(loc, "outs/", dataset, "/", dataset, "_", "rand_2.pdf"), plot = grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2), height = 8, width = 12)

}