### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("dplyr", "foreach", "doParallel", "MODIS", "raster", "rgdal", 
         "hydroGOF", "lubridate", "caret")
Orcs::loadPkgs(lib)

## load functions
source("R/movingModel.R")
source("R/funs.R")

## parallelization
cl1 <- makeCluster(detectCores() / 2)
registerDoParallel(cl1)


### lst -----

## import files (daytime)
fls_lst <- list.files("data/MOD11A2.005/qc", 
                      pattern = "^MOD11A2.*_Day_1km.tif$", full.names = TRUE)

fls_lst <- fls_lst[grep("2011", fls_lst)[1]:grep("2015361", fls_lst)]

rst_lst <- stack(fls_lst) - 273.15
dts_lst <- MODIS::extractDate(fls_lst)$inputLayerDates

# ## import files (nighttime)
# fls_lst_night <- list.files("data/MOD11A2.005/gf", full.names = TRUE,
#                             pattern = "^MOD11A2.*_Night_1km.tif$")
# 
# fls_lst_night <- fls_lst_night[grep("2010", fls_lst_night)[1]:length(fls_lst_night)]
# 
# rst_lst_night <- stack(fls_lst_night) - 273.15
# dts_lst_night <- MODIS::extractDate(fls_lst_night)$inputLayerDates


### ndvi -----

## import files
fls_ndvi <- list.files("data/MCD09Q1.006/ndvi", 
                       pattern = "^MCD09Q1.*.tif$", full.names = TRUE)

fls_ndvi <- fls_ndvi[grep("2011", fls_ndvi)[1]:grep("2015361", fls_ndvi)]

rst_ndvi <- stack(fls_ndvi)
dts_ndvi <- MODIS::extractDate(fls_ndvi, 24, 30)$inputLayerDates

# ## remove unavailable lst files
# rst_lst <- rst_lst[[which(dts_lst %in% dts_ndvi)]]
# rst_lst_night <- rst_lst_night[[which(dts_lst_night %in% dts_ndvi)]]

## resample ndvi
rst_ndvi_res <- resample(rst_ndvi, rst_lst)


### temperature-vegetation index (tvx) related slope, intercept, and
### regression coefficient -----

## target folder
dir_tvx <- "data/MOD11A2.005/tvx_7by7"
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
    
    mod <- movingModel(rst_lst[[i]], rst_ndvi_res[[i]], 
                       directions = movingWindow(7L))
    mod <- stack(mod)
    
    # save metrics
    writeRaster(mod, filename = fls_tvx, format = "GTiff", overwrite = TRUE)
  }
}

## convert to matrices
rst_intercept <- stack(lapply(lst_tvx, "[[", 2))
mat_intercept <- as.matrix(rst_intercept)

rst_slope <- stack(lapply(lst_tvx, "[[", 1))
mat_slope <- as.matrix(rst_slope)


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


### evaluation -----

## loop over plots
out <- foreach(h = as.character(unique(dat_plots$plotID)), 
               .packages = lib, .combine = "rbind") %dopar% {

  ## temporary parallel cluster               
  cl2 <- makeCluster(2)
  registerDoParallel(cl2)
                 
  ## extract daily maximum temperature
  dat_plots %>%
    filter(plotID == h & !is.na(Ta_200)) %>%
    group_by(Date = as.Date(datetime)) %>%
    filter(length(Ta_200) == 24) %>%
    summarise(Ta_200 = max(Ta_200, na.rm = TRUE)) %>%
    data.frame() -> dat_sub
  
  ## synchronise with modis time steps, i.e. aggregate 8-day mean maximum values 
  dat_agg <- foreach(i = unique(lubridate::year(dat_sub$Date)), 
                     .combine = "rbind") %do% {

    # merge available dates with continuous daily time series
    dts0 <- seq(as.Date(paste0(i, "-01-01")), as.Date(paste0(i, "-12-31")), 1)
    mrg <- merge(data.frame(Date = dts0), dat_sub, all.x = TRUE)
    
    # aggregate 8-day intervals
    mrg$Cuts <- cut(mrg$Date, "8 days")
    mrg <- mrg[complete.cases(mrg), ]
    
    mrg %>% 
      group_by(Cuts) %>%
      summarise(Ta_200 = round(mean(Ta_200, na.rm = TRUE), 2)) %>%
      data.frame()
  }
  
  ## extract intercept and slope values
  dts_lst <- MODIS::extractDate(names(rst_lst))$inputLayerDates
  avl <- dts_lst %in% strftime(dat_agg$Cuts, "%Y%j")
  
  if (all(!avl)) {
    return(data.frame(PlotID = h, NDVImax = NA, rsq = NA, rmse = NA, n = NA))
  }
  
  metrics <- foreach(i = 2:1, .export = ls(envir = environment()), 
                     .packages = lib) %dopar% {
    rst <- raster::stack(sapply(lst_tvx, "[[", i))
    rst <- rst[[which(avl)]]
    
    val <- as.numeric(raster::extract(rst, plots[plots@data$PlotID == h, ]))
    
    if (all(is.na(val))) {
      data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates, 
                 rep(NA, length(MODIS::extractDate(names(rst))$inputLayerDates)))
    } else {
      data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates, val)
    }
  }
  
  metrics <- Reduce(function(...) merge(..., by = "Date", all = TRUE), metrics)
  names(metrics)[2:3] <- c("Intercept", "Slope")
  
  if (all(is.na(metrics$Intercept)) | all(is.na(metrics$Slope))) {
    return(data.frame(PlotID = h, NDVImax = NA, rsq = NA, rmse = NA, n = NA))
  }
  
  ## loop over sequence of ndvi_max values
  dat_agg$Date <- strftime(dat_agg$Cuts, "%Y%j")
  dat_mrg <- merge(dat_agg, metrics, all.y = TRUE)
  
  n <- sum(complete.cases(dat_mrg))
  
  ta_obs <- dat_mrg$Ta_200
  
  stats <- foreach(i = seq(.2, 1, .01), .combine = "rbind") %do% {
    
    # estimate temperature
    ta_prd <- (dat_mrg$Slope * i + dat_mrg$Intercept)
    
    # calculate r-squared, rmse, and bias
    data.frame(PlotID = h, NDVImax = i, 
               rsq = summary(lm(ta_prd ~ ta_obs))$r.squared, 
               rmse = rmse(ta_obs, ta_prd), 
               bias = pbias(ta_prd, ta_obs), n = n)
  }
  
  parallel::stopCluster(cl2)
  stats[which.min(stats$rmse), ]
}

saveRDS(out, file = "data/results/ndvimax_rmse_7by7_terra.rds")
