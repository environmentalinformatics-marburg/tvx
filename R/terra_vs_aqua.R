out_terra <- readRDS("data/results/ndvimax_rmse_(terra).rds")
out_aqua <- readRDS("data/results/ndvimax_rmse.rds")

out_terra[[1]]
out_aqua[[1]]

id <- sapply(sapply(out_terra, "[[", 1), function(i) class(i) == "data.frame")
lst_terra <- lapply(out_terra, "[[", 1)[id]
dat_terra <- do.call("rbind", lst_terra)

id <- sapply(sapply(out_aqua, "[[", 1), function(i) class(i) == "data.frame")
lst_aqua <- lapply(out_aqua, "[[", 1)[id]
dat_aqua <- do.call("rbind", lst_aqua)


xyplot(dat_terra$rmse ~ dat_aqua$rmse)
