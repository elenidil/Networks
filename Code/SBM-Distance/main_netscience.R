rm(list=ls(all=TRUE))

loc     <- "/soe/jsosamar/NETS/"
dataset <- "netscience"

K_range <- c(2, 4, 8)
Q_range <- c(2, 4, 8)
n_sams  <- 75000
n_burn  <- 25000
n_skip  <- 1

path_outs <- paste0(loc, "outs/", dataset, "/") 
if (!dir.exists(path_outs)) dir.create(path_outs)

source(paste0(loc, "code/wrapper.R"))

wrapper_0(                  n_sams, n_burn, n_skip, model = "erdos"  , dataset, loc)
wrapper_1(K_range,          n_sams, n_burn, n_skip, model = "dist"   , dataset, loc)
wrapper_1(K_range,          n_sams, n_burn, n_skip, model = "class"  , dataset, loc)
wrapper_1(K_range,          n_sams, n_burn, n_skip, model = "eigen"  , dataset, loc)
wrapper_2(K_range, Q_range, n_sams, n_burn, n_skip, model = "clist"  , dataset, loc)
wrapper_2(K_range, Q_range, n_sams, n_burn, n_skip, model = "cleigen", dataset, loc)
wrapper_2(K_range, Q_range, n_sams, n_burn, n_skip, model = "class2" , dataset, loc)

get_results(dataset, loc)
plot_measures(K_range, Q_range, dataset, loc)