### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("dplyr", "foreach", "doParallel", "MODIS", "raster", "rgdal", "grid",
         "hydroGOF", "lubridate", "caret", "latticeExtra", "reshape2", "Rsenal")
Orcs::loadPkgs(lib)

## load functions
source("R/movingModel.R")
source("R/funs.R")

## parallelization
cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)


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


### temperature-vegetation index (tvx) related slope, intercept, and
### regression coefficient -----

## dem
rst_dem <- raster("data/dem/DEM_ARC1960_30m_Hemp.tif")
rst_dem_res <- resample(rst_dem, rst_ndvi_res)

## target folder
dir_tvx <- "data/MOD11A2.005/tvx_7by7_enhanced"
if (!dir.exists(dir_tvx)) dir.create(dir_tvx)

## loop over layers
lst_tvx <- foreach(i = 1:nlayers(rst_lst), .packages = "raster") %dopar% {
  
  # target file
  fls_tvx <- paste0(dir_tvx, "/", gsub("LST", "TVX", names(rst_lst)[i]), ".tif")
  
  # import file if available
  if (file.exists(fls_tvx)) {
    stack(fls_tvx)
    
    # else fit moving window model  
  } else {
    mod <- movingModel(rst_lst[[i]], rst_ndvi_res[[i]], rst_dem_res,
                       directions = movingWindow(7L))
    mod <- stack(mod)
    
    # save metrics
    writeRaster(mod, filename = fls_tvx, format = "GTiff", overwrite = TRUE)
  }
}

## convert to matrices
rst_intercept <- stack(lapply(lst_tvx, "[[", 3))
mat_intercept <- as.matrix(rst_intercept)

rst_slope1 <- stack(lapply(lst_tvx, "[[", 1)) # lst
mat_slope1 <- as.matrix(rst_slope1)
rst_slope2 <- stack(lapply(lst_tvx, "[[", 2)) # dem
mat_slope2 <- as.matrix(rst_slope2)


### station data -----

## import plot coordinates
plots <- readOGR(dsn = "data/station_data", 
                 layer = "PlotPoles_ARC1960_mod_20140807_final", 
                 p4s = "+init=epsg:21037")

plots <- subset(plots, PoleType == "AMP")

## import available plot data 
fls_plots <- list.files("data/station_data/plots", pattern = ".csv", 
                        full.names = TRUE)  

dat_plots <- foreach(i = fls_plots, .combine = "rbind") %dopar% {
  read.csv(i)[, c("plotID", "datetime", "Ta_200")]
}

dat_plots$habitat <- substr(dat_plots$plotID, 1, 3)
dat_plots$habitat[dat_plots$plotID == "mch0"] <- "fpo"
dat_plots$habitat[dat_plots$plotID == "mwh0"] <- "fpd"


### select training plots -----

## loop over plots
dat_agg <- foreach(h = as.character(unique(dat_plots$plotID)),
                      .packages = lib, .combine = "rbind") %dopar% {
                        
  ## extract daily maximum temperature
  dat_plots %>%
    filter(plotID == h & !is.na(Ta_200)) %>%
    group_by(Date = as.Date(datetime)) %>%
    filter(length(Ta_200) == 24) %>%
    summarise(Ta_200 = max(Ta_200, na.rm = TRUE)) %>%
    data.frame() -> dat_sub
  
  data.frame(PlotID = h, Ta_200 = round(mean(dat_sub$Ta_200, na.rm = TRUE), 2))                                            
}

dat_mrg <- merge(plots@data, dat_agg, by = "PlotID", all.y = TRUE)


clr <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))

p_below <- xyplot(Z_DEM_HMP ~ Ta_200, data = subset(dat_mrg, Z_DEM_HMP <= 2320), 
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, col = "grey75", ...)
                    panel.ablineq(lm(y ~ x), r.squared = TRUE, rotate = TRUE, 
                                  lty = 2, lwd = 2)
                  }, xlim = c(8, 36), ylim = c(750, 4750))

p_above <- xyplot(Z_DEM_HMP ~ Ta_200, data = subset(dat_mrg, Z_DEM_HMP > 2320), 
                  panel = function(x, y, ...) {
                    panel.xyplot(x, y, col = "grey75", ...)
                    panel.ablineq(lm(y ~ x), r.squared = TRUE, rotate = TRUE, 
                                  lty = 2, lwd = 2)
                  })

p_below + 
  as.layer(p_above)

