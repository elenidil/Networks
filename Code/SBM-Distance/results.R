rm(list=ls(all=TRUE))

loc     <- "C:/Users/Juan Camilo/Dropbox/PAPERS/NETS/"
dataset <- "lesmis"

K_range <- c(2, 4, 8)
Q_range <- c(2, 4, 8)
n_sams  <- 75000

path_outs <- paste0(loc, "outs/", dataset, "/")
if (!dir.exists(path_outs)) dir.create(path_outs)

source(paste0(loc, "code/wrapper.R"))

plot_measures2(K_range, Q_range, dataset, loc)
plot_stats2(K_range, Q_range, dataset, loc)
plot_rand2(K_range, Q_range, dataset, loc, n_sams)

#zach
#dol
#fblog
#lazega
#indus
#jazz
#lesmis