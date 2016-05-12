## packages
library(lattice)
library(RColorBrewer)
library(dplyr)

## import data
dat_terra_5by5 <- readRDS("data/results/ndvimax_rmse_5by5_terra.rds")
dat_terra_7by7 <- readRDS("data/results/ndvimax_rmse_7by7_terra.rds")

## concat datasets
getMetrics <- function(x) {
  lst <- lapply(x, "[[", 1)
  
  id <- sapply(lst, class) == "data.frame"
  lst <- lst[id]
  
  do.call("rbind", lst)
}

lst_terra <- lapply(list(dat_terra_5by5, dat_terra_7by7), getMetrics)

for (i in 1:2) {
  lst_terra[[i]]$Group <- ifelse(i == 1, "5 x 5", "7 x 7")
}

dat_terra <- do.call("rbind", lst_terra)
dat_terra <- reshape2::melt(dat_terra, id.vars = c("PlotID", "Group"))

## display density distributions
clr <- brewer.pal(4, "PuOr")[c(1, 4)]
names(clr) <- c("5 x 5", "7 x 7")

lty <- c("5 x 5" = 1, "7 x 7" = 2)

densityplot( ~ value | variable, data = dat_terra, groups = Group, lty = lty,
             scales = list(x = list(relation = "free")), col = clr, lwd = 2,
             layout = c(3, 1), panel = function(x, ...) {
               panel.densityplot(x, ..., plot.points = FALSE)
             }, key = list("lines" = list(col = clr, lty = lty), 
                           "text" = list(c("5 x 5", "7 x 7")), columns = 2L))

## summarise metrics (ndvi_max, rsq, rmse) per matrix group
dat_terra %>% 
  group_by(Group, variable) %>%
  summarise(mean(value, na.rm = TRUE))


cof1, 