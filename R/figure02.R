### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("latticeExtra", "reshape2", "grid")
Orcs::loadPkgs(lib)


### plot prediction -----

## import data
dat_training <- readRDS("data/results/stats_intra_7by7_terra.rds")

## reorder factor levels
dat_training$habitat <- factor(dat_training$habitat, 
                               levels = rev(sortElevation(df = FALSE)))

## melt data
prm1 <- c("habitat", "TrainRsq", "TrainRMSE", "TrainRMSEse", 
         "Rsq", "RMSE", "RMSEse")
training_stats_mlt <- melt(dat_training[, prm1])

## strip levels
lvl1 <- c(expression("a)" ~ bold("R-squared")), expression("b)" ~ bold("RMSE")))

## colors
clr <- RColorBrewer::brewer.pal(3, "PuOr")

## axis limits
xlim_rsq <- c(0, .65); xlim_rmse <- c(0, 12.5)

## create plot
p1 <- dotplot(habitat ~ value | variable, xlim = list(xlim_rsq, xlim_rmse), 
              data = subset(training_stats_mlt, variable %in% c("Rsq", "RMSE")),
              cex = 1.2, xlab = NULL, ylab = NULL, 
              strip = strip.custom(factor.levels = lvl1, bg = clr[1]),
              scales = list(x = list(relation = "free", labels = NULL)), 
              par.settings = list(layout.heights = list(strip = 1.1)),
              panel = function(x, y, ...) {
                panel.dotplot(x, y, col.line = "grey90")
                for (i in 1:length(levels(dat_training$habitat))) {
                  hab <- levels(dat_training$habitat)[i]
                  tmp <- subset(dat_training, habitat == hab)
                  
                  if (panel.number() == 2) {
                    xmn <- tmp$RMSE - tmp$RMSEse
                    xmn <- ifelse(xmn < 0, 0, xmn)
                    
                    xmx <- tmp$RMSE + tmp$RMSEse
                    
                    panel.lines(x = c(xmn, xmx), y = i, col = "black", 
                                lwd = 1.2)
                    panel.dotplot(x = c(xmn, xmx), y = i, pch = "|", 
                                  cex = 1.1, col.line = "transparent", 
                                  col = "black")
                  }
                }
                
                panel.dotplot(x, y, col = "black", col.line = "transparent", 
                              ...)
              })


### habitat prediction -----

## import data
prediction_stats <- readRDS("data/results/stats_inter_7by7_terra.rds")

## reorder factor levels
prediction_stats$habitat <- factor(prediction_stats$habitat, 
                                   levels = rev(sortElevation(df = FALSE)))

## melt data
prm2 <- c("habitat", "rsq", "rmse", "rmse.se")
prediction_stats_mlt <- melt(prediction_stats[, prm2])

## strip levels
lvl2 <- c(expression("c)" ~ bold("R-squared")), expression("d)" ~ bold("RMSE")))

## create plot
p2 <- dotplot(habitat ~ value | variable, xlim = list(xlim_rsq, xlim_rmse), 
              data = subset(prediction_stats_mlt, variable %in% c("rsq", "rmse")),
              cex = 1.2, col = "black", xlab = NULL, ylab = NULL, 
              strip = strip.custom(factor.levels = lvl2, bg = clr[3]),
              scales = list(x = list(relation = "free", 
                                     at = list(seq(0, .6, .1), seq(0, 12, 2)), 
                                     labels = list(c("0", "", "0.2", "", "0.4", "", "0.6"), 
                                                   c("0", "", "4", "", "6", "", "8", "", "12")))), 
              par.settings = list(layout.heights = list(strip = 1.1)),
              panel = function(x, y, ...) {
                panel.dotplot(x, y, col.line = "grey90")
                for (i in 1:length(levels(prediction_stats$habitat))) {
                  hab <- levels(prediction_stats$habitat)[i]
                  tmp <- subset(prediction_stats, habitat == hab)
                  
                  if (panel.number() == 2) {
                    xmn <- tmp$rmse - tmp$rmse.se
                    xmn <- ifelse(xmn < 0, 0, xmn)
                    
                    xmx <- tmp$rmse + tmp$rmse.se
                    
                    panel.lines(x = c(xmn, xmx), y = i, col = "black", 
                                lwd = 1.2)
                    panel.dotplot(x = c(xmn, xmx), y = i, pch = "|", cex = 1.1,
                                  col = "black", col.line = "transparent")
                  }
                }
                
                panel.dotplot(x, y, col.line = "transparent", ...)
              })


### visualize -----

## standalone .tiff version
tiff("vis/figure02.tiff", width = 19, height = 16, units = "cm", res = 600, 
     compression = "lzw")
grid.newpage()

## add training stats
vp0 <- viewport(x = .05, y = 1, width = .95, height = .5, just = c("left", "top"))
pushViewport(vp0)
print(p1, newpage = FALSE)

## add test stats
upViewport()
vp1 <- viewport(x = .05, y = .636, width = .95, height = .5, just = c("left", "top"))
pushViewport(vp1)
print(p2, newpage = FALSE)

## add y-axis label
upViewport()
grid.text("Habitat type", x = .025, y = .58, rot = 90)

dev.off()
