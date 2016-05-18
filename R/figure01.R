### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("foreach", "MODIS", "Rsenal", "latticeExtra", "reshape2", "grid", 
         "gridBase", "mapview")
Orcs::loadPkgs(lib)

## load functions
source("R/movingModel.R")
source("R/funs.R")
source("R/visKili.R")


### lst -----

## import files (daytime)
fls_lst <- list.files("data/MOD11A2.005/qc", 
                      pattern = "^MOD11A2.*_Day_1km.tif$", full.names = TRUE)

fls_lst <- fls_lst[grep("2011001", fls_lst)[1]:grep("2016121", fls_lst)]

rst_lst <- stack(fls_lst) - 273.15
mat_lst <- as.matrix(rst_lst)
dts_lst <- MODIS::extractDate(fls_lst)$inputLayerDates


### ndvi -----

## import files
fls_ndvi <- list.files("data/MCD09Q1.006/ndvi", 
                       pattern = "^MCD09Q1.*.tif$", full.names = TRUE)

fls_ndvi <- fls_ndvi[grep("2011001", fls_ndvi)[1]:grep("2016121", fls_ndvi)]

rst_ndvi <- stack(fls_ndvi)
dts_ndvi <- MODIS::extractDate(fls_ndvi)$inputLayerDates

## remove unavailable lst files
rst_ndvi <- rst_ndvi[[-which(!dts_ndvi %in% dts_lst)]]

## resample ndvi
rst_ndvi_res <- resample(rst_ndvi, rst_lst)
mat_ndvi_res <- raster::as.matrix(rst_ndvi_res)


### study area -----

# rst_kili <- kiliAerial(projection = "+init=epsg:4326", type = "google")
# rst_kili <- writeRaster(rst_kili, "data/kili_ll", format = "GTiff", 
#                         overwrite = TRUE)

rst_kili <- stack("data/kili_ll.tif")


### correlation -----

# loop over matrix rows (i.e., pixels)
mat_cor <- foreach(k = 1:nrow(mat_lst), .combine = "rbind") %do% {
  
  # if enough valid values are available, compute linear regression metrics
  dat <- data.frame(y = mat_lst[k, ], x = mat_ndvi_res[k, ])
  
  if (sum(complete.cases(dat)) >= 23) {
    dat <- dat[complete.cases(dat), ]
    
    r <- cor(dat$x, dat$y)
    mod <- lm(y ~ x, data = dat)
    
    # return metrics only if highly significant (p < 0.001)
    data.frame(n = sum(complete.cases(dat)), 
               r = r, p = Orcs::pvalue(mod))
    
  } else {
    return(NA)
  }
}

# write r values to raster
rst_r <- trim(projectRaster(setValues(rst_lst[[1]], mat_cor[, 2]), 
                            crs = "+init=epsg:4326"))
rst_p <- trim(projectRaster(setValues(rst_lst[[1]], mat_cor[, 3]), 
                            crs = "+init=epsg:4326"))

rst_ndvi_mu <- trim(projectRaster(calc(rst_ndvi_res, fun = function(x) {
  mean(x, na.rm = TRUE)
}), crs = "+init=epsg:4326"))

rst_robust <- overlay(rst_r, rst_ndvi_mu, fun = function(x, y) {
  x[y[] < 0.15] <- NA
  return(x)
})

rst_r_sig <- overlay(rst_robust, rst_p, fun = function(x, y) {
  x[y[] >= .001] <- NA
  return(x)
})

clr <- Rsenal::envinmrPalette(200)

cnt_r_sig <- rst_r_sig

# ## smooth significant areas, discard majority of island pixels
# cnt_r_sig <- focal(cnt_r_sig, fun = function(x) {
#   if (sum(is.na(x)) > 3) NA else mean(x, na.rm = TRUE)
# }, w = matrix(1/9, nc = 3, nr = 3), pad = FALSE)
# cnt_r_sig[1, ] <- rst_r_sig[1, ]
# cnt_r_sig[, 1] <- rst_r_sig[, 1]
# cnt_r_sig[nrow(cnt_r_sig), ] <- rst_r_sig[nrow(rst_r_sig), ]
# cnt_r_sig[, ncol(cnt_r_sig)] <- rst_r_sig[, ncol(rst_r_sig)]

cnt_r_sig[is.na(cnt_r_sig[])] <- -999
cnt_r_sig[cnt_r_sig[] > -999] <- NA
cnt_r_sig <- overlay(cnt_r_sig, rst_ndvi_mu, fun = function(x, y) {
  x[y[] < 0.15 | is.na(y[])] <- NA
  return(x)
})

spy_r_sig <- rasterToPolygons(cnt_r_sig)
spy_r_dis <- rasterToPolygons(cnt_r_sig, dissolve = TRUE)


### visualize -----

## study area
scale <- list("SpatialPolygonsRescale", layout.scale.bar(), scale = 0.08998623, 
              offset = c(37.05, -3.38), fill = c("transparent", "black"))
text1 = list("sp.text", c(37.05, -3.36), "0", cex = .5, font = 2)
text2 = list("sp.text", c(37.16, -3.36), "10 km", cex = .5, font = 2)

arrow <- list("SpatialPolygonsRescale", layout.north.arrow(type = 1), 
              offset = c(37, -3.41), scale = .075)

p_kili <- spplot(rst_robust, col.regions = "transparent", 
                 scales = list(draw = TRUE, cex = .6), colorkey = FALSE, 
                 sp.layout = list(rgb2spLayout(rst_kili, c(.01, .998)), 
                                  scale, text1, text2, 
                                  list("sp.text", loc = c(37.02, -2.86), 
                                       txt = "a)", font = 2, cex = .6, 
                                       adj = c(.1, 1), col = "black")))

## topographic map
p_topo <- visKili(cex = .6, lwd = .1, ext = rst_kili)

## standalone .tif version
tiff("vis/figure01.tiff", width = 19, height = 9, units = "cm", res = 300, 
     compression = "lzw")
plot.new()

# add image of study area
vp0 <- viewport(-.025, 0, .55, 1, just = c("left", "bottom"))
pushViewport(vp0)
print(p_kili, newpage = FALSE)

# add topographic map
downViewport(trellis.vpname("figure"))
vp_topo <- viewport(x = .65, y = .6, just = c("left", "bottom"), 
                    width = .425, height = .475)
pushViewport(vp_topo)
print(p_topo, newpage = FALSE)

# add raster with r values
upViewport(0)
vp1 <- viewport(.425, 0, .55, 1, just = c("left", "bottom"))
pushViewport(vp1)
print(spplot(rst_robust, 
             at = seq(-.85, .85, .01), col.regions = clr, colorkey = FALSE, 
             scales = list(draw = TRUE, cex = .6, y = list(col = "transparent")), 
             par.settings = list(panel.background=list(col="grey80")), 
             sp.layout = list("sp.text", loc = c(37.02, -2.86), 
                              txt = "b)", font = 2, cex = .6, 
                              adj = c(.1, 1), col = "black")), 
      newpage = FALSE)

# add colorkey
downViewport(trellis.vpname("figure"))

vp2 <- viewport(1.075, 0, .1, 1, just = c("left", "bottom"))
pushViewport(vp2)
draw.colorkey(key = list(labels = list(cex = .6), col = clr, 
                         at = seq(-.85, .85, .01), width = .6, height = .75, 
                         space = "right"), draw = TRUE)

grid.text("r", 0, .95, gp = gpar(fontface = "bold", cex = .7))

# add hatches
upViewport()
par(new = TRUE, fig = gridFIG(), mai = c(0, 0, 0, 0))
plot(spy_r_sig, border = "transparent", bg = "transparent", density = 10, 
     col = "grey50", xlim = c(xmin(rst_r_sig), xmax(rst_r_sig)), xaxs = "i", 
     yaxs = "i", ylim = c(ymin(rst_r_sig), ymax(rst_r_sig)))

# add contour lines
par(new = TRUE, fig = gridFIG(), mai = c(0, 0, 0, 0))
plot(spy_r_dis, border = "grey30", bg = "transparent", lwd = 1.2,
     xlim = c(xmin(rst_r_sig), xmax(rst_r_sig)), xaxs = "i", 
     yaxs = "i", ylim = c(ymin(rst_r_sig), ymax(rst_r_sig)))

# add black margin
grid.rect(gp = gpar(fill = "transparent"))

dev.off()
