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
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### lst -----

## import files (daytime)
fls_lst <- list.files("data/MOD11Q1.006", 
                      pattern = "Day_1km.tif$", full.names = TRUE)
rst_lst <- stack(fls_lst)
mat_lst <- as.matrix(rst_lst)
dts_lst <- MODIS::extractDate(fls_lst)$inputLayerDates


### ndvi -----

## list available files
fls_ndvi <- list.files("data/MCD13Q1.006/whittaker",
                       pattern = "NDVI.tif$", full.names = TRUE)
dts_ndvi <- MODIS::extractDate(fls_ndvi)$inputLayerDates

## remove unavailable lst files
fls_ndvi <- fls_ndvi[-which(!dts_ndvi %in% dts_lst)]
rst_ndvi <- stack(fls_ndvi)
mat_ndvi <- as.matrix(rst_ndvi)


### temperature-vegetation index (tvx) related slope, intercept, and
### regression coefficient -----

## target folder
dir_tvx <- "data/MOD11Q1.006/tvx_7by7"
if (!dir.exists(dir_tvx)) dir.create(dir_tvx)

## loop over layers
fls_tvx <- paste0(dir_tvx, "/", basename(fls_lst))
fls_tvx <- gsub("LST", "TVX", fls_tvx)

lst_tvx <- foreach(i = 1:nlayers(rst_lst), .packages = "raster") %dopar% {
  
  # import file if available
  if (file.exists(fls_tvx[i])) {
    stack(fls_tvx[i])
    
    # else fit moving window model  
  } else {
    
    mod <- movingModel(rst_lst[[i]], rst_ndvi[[i]], 
                       directions = movingWindow(7L))
    mod <- stack(mod)
    
    # save metrics
    writeRaster(mod, filename = fls_tvx[i], format = "GTiff", overwrite = TRUE)
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

# ## remove plots in immediate vicinity to water (chala) or settlements (moshi)
# rst_lc <- raster("data/MCD12Q1.051/MCD12Q1.A2013001.Land_Cover_Type_5.tif")
# 
# id <- extract(rst_lc, plots, cellnumbers = TRUE)[, 1]
# 
# invalid <- foreach(i = id, .combine = "c", .packages = "raster") %dopar% {
#   cells <- adjacent(rst_lc, i, directions = movingWindow(15L), pairs = FALSE)
#   any(rst_lc[cells] %in% c(0, 9))
# }

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

# ## loop over plots
# dat_select <- foreach(h = as.character(unique(dat_plots$plotID)),
#                .packages = lib, .combine = "rbind") %dopar% {
# 
#   ## extract daily maximum temperature
#   dat_plots %>%
#     filter(plotID == h & !is.na(Ta_200)) %>%
#     group_by(Date = as.Date(datetime)) %>%
#     filter(length(Ta_200) == 24) %>%
#     summarise(Ta_200 = max(Ta_200, na.rm = TRUE)) %>%
#     data.frame() -> dat_sub
# 
#   ## synchronise with modis time steps, i.e. aggregate 8-day mean maximum values
#   dat_agg <- foreach(i = unique(lubridate::year(dat_sub$Date)),
#                      .combine = "rbind") %do% {
# 
#     # merge available dates with continuous daily time series
#     dts0 <- seq(as.Date(paste0(i, "-01-01")), as.Date(paste0(i, "-12-31")), 1)
#     mrg <- merge(data.frame(Date = dts0), dat_sub, all.x = TRUE)
# 
#     # aggregate 8-day intervals
#     mrg$Cuts <- cut(mrg$Date, "8 days")
#     mrg <- mrg[complete.cases(mrg), ]
# 
#     mrg %>%
#       group_by(Cuts) %>%
#       summarise(Ta_200 = round(mean(Ta_200, na.rm = TRUE), 2)) %>%
#       data.frame()
#   }
# 
#   ## extract intercept and slope values
#   dts_lst <- MODIS::extractDate(names(rst_lst))$inputLayerDates
#   avl <- dts_lst %in% strftime(dat_agg$Cuts, "%Y%j")
# 
#   if (all(!avl) | sum(avl) == 1) {
#     return(data.frame(PlotID = h, NDVImax = NA, rsq = NA, rmse = NA, bias = NA,
#                       n = ifelse(sum(avl) == 1, 1, NA)))
#   }
# 
#   metrics <- foreach(i = list(mat_intercept, mat_slope),
#                      j = list(rst_intercept, rst_slope)) %do% {
# 
#     mat <- i[, avl]
#     rst <- j[[which(avl)]]
# 
#     cell_id <- cellFromXY(j, plots[plots@data$PlotID == h, ])
#     val <- mat[cell_id, ]
# 
#     if (all(is.na(val))) {
#       data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates,
#                  rep(NA, length(MODIS::extractDate(names(rst))$inputLayerDates)))
#     } else {
#       data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates, val)
#     }
#   }
# 
#   metrics <- Reduce(function(...) merge(..., by = "Date", all = TRUE), metrics)
#   names(metrics)[2:3] <- c("Intercept", "Slope")
# 
#   if (all(is.na(metrics$Intercept)) | all(is.na(metrics$Slope))) {
#     return(data.frame(PlotID = h, NDVImax = NA, rsq = NA, rmse = NA, bias = NA, n = NA))
#   }
# 
#   ## loop over sequence of ndvi_max values
#   dat_agg$Date <- strftime(dat_agg$Cuts, "%Y%j")
#   dat_mrg <- merge(dat_agg, metrics, all.y = TRUE)
# 
#   n <- sum(complete.cases(dat_mrg))
# 
#   ta_obs <- dat_mrg$Ta_200
# 
#   stats <- foreach(i = seq(.2, 1, .01), .combine = "rbind") %do% {
# 
#     # estimate temperature
#     ta_prd <- (dat_mrg$Slope * i + dat_mrg$Intercept)
# 
#     # calculate r-squared, rmse, and bias
#     data.frame(PlotID = h, NDVImax = i,
#                rsq = summary(lm(ta_prd ~ ta_obs))$r.squared,
#                rmse = rmse(ta_obs, ta_prd),
#                bias = pbias(ta_prd, ta_obs), n = n)
#   }
# 
#   # parallel::stopCluster(cl2)
#   stats[which.min(stats$rmse), ]
# }
# 
# dat_select <- dat_select[complete.cases(dat_select), ]
# 
# saveRDS(dat_select, file = "data/results/ndvimax_rmse_7by7_area_terra.rds")
dat_select <- readRDS("data/results/ndvimax_rmse_7by7_area_terra.rds")

## extract habitat information
dat_select$habitat <- substr(dat_select$PlotID, 1, 3)
dat_select$habitat[dat_select$PlotID == "mch0"] <- "fpo"
dat_select$habitat[dat_select$PlotID == "mwh0"] <- "fpd"

dat_select %>%
  filter(n >= 10) %>%
  mutate(ratio = rmse / n * (-1)) %>%
  group_by(habitat) %>%
  top_n(n = 2, wt = ratio) %>%
  data.frame() -> dat_sel


### evaluate approach at selected training plots -----

# ## loop over selected training plots
# dat_training <- foreach(h = unique(dat_sel$habitat),
#                         .packages = lib, .combine = "rbind") %dopar% {
# 
#   dat_sel_sub <- subset(dat_sel, habitat == h)
# 
#   dat_mrg <- foreach(k = as.character(dat_sel_sub$PlotID),
#                      .combine = "rbind") %do% {
# 
#     ## extract daily maximum temperature
#     dat_plots %>%
#       filter(plotID == k & !is.na(Ta_200)) %>%
#       group_by(Date = as.Date(datetime)) %>%
#       filter(length(Ta_200) == 24) %>%
#       summarise(Ta_200 = max(Ta_200, na.rm = TRUE)) %>%
#       data.frame() -> dat_sub
# 
#     ## synchronise with modis time steps, i.e. aggregate 8-day mean maximum values
#     dat_agg <- foreach(i = unique(lubridate::year(dat_sub$Date)),
#                        .combine = "rbind") %do% {
# 
#       # merge available dates with continuous daily time series
#       dts0 <- seq(as.Date(paste0(i, "-01-01")), as.Date(paste0(i, "-12-31")), 1)
#       mrg <- merge(data.frame(Date = dts0), dat_sub, all.x = TRUE)
# 
#       # aggregate 8-day intervals
#       mrg$Cuts <- cut(mrg$Date, "8 days")
#       mrg <- mrg[complete.cases(mrg), ]
# 
#       mrg %>%
#         group_by(Cuts) %>%
#         summarise(Ta_200 = round(mean(Ta_200, na.rm = TRUE), 2)) %>%
#         data.frame()
#     }
# 
#     ## extract intercept and slope values
#     dts_lst <- MODIS::extractDate(names(rst_lst))$inputLayerDates
#     avl <- dts_lst %in% strftime(dat_agg$Cuts, "%Y%j")
# 
#     if (all(!avl)) {
#       return(data.frame(PlotID = k, NDVImax = NA, rsq = NA, rmse = NA))
#     }
# 
#     metrics <- foreach(i = list(mat_intercept, mat_slope),
#                        j = list(rst_intercept, rst_slope)) %do% {
# 
#                          mat <- i[, avl]
#                          rst <- j[[which(avl)]]
# 
#                          cell_id <- cellFromXY(j, plots[plots@data$PlotID == k, ])
#                          val <- mat[cell_id, ]
# 
#                          if (all(is.na(val))) {
#                            data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates,
#                                       rep(NA, length(MODIS::extractDate(names(rst))$inputLayerDates)))
#                          } else {
#                            data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates, val)
#                          }
#                        }
# 
#     metrics <- Reduce(function(...) merge(..., by = "Date", all = TRUE), metrics)
#     names(metrics)[2:3] <- c("Intercept", "Slope")
# 
#     if (all(is.na(metrics$Intercept)) | all(is.na(metrics$Slope))) {
#       return(data.frame(PlotID = k, NDVImax = NA, rsq = NA, rmse = NA))
#     }
# 
#     ## loop over sequence of ndvi_max values
#     dat_agg$Date <- strftime(dat_agg$Cuts, "%Y%j")
#     merge(dat_agg, metrics, all.y = TRUE)
#   }
# 
#   n <- sum(complete.cases(dat_mrg))
# 
#   ta_obs <- dat_mrg$Ta_200
# 
#   tst <- data.frame(Ta_obs = ta_obs, Slope = dat_mrg$Slope, Intercept = dat_mrg$Intercept)
#   tst <- tst[complete.cases(tst), ]
# 
#   split <- .5
# 
#   res <- do.call(function(...) colMeans(rbind(...), na.rm = TRUE),
#                  lapply(1:100, function(g) {
# 
#     set.seed(g)
#     trainIndex <- sort(sample(1:length(tst$Ta_obs),
#                               floor(length(tst$Ta_obs) * split)))
# 
#     data_train <- tst[trainIndex, ]
#     data_test <- tst[-trainIndex, ]
# 
#     ## determine ndvi_max
#     stats <- do.call("rbind", lapply(seq(.2, 1, .01), function(i) {
# 
#       # estimate temperature
#       ta_prd <- (data_train$Slope * i + data_train$Intercept)
# 
#       stats_trn <- Rsenal::regressionStats(ta_prd, data_train$Ta_obs, adj.rsq = FALSE)
# 
#       # calculate r-squared, rmse
#       data.frame(PlotID = h, NDVImax = i,
#                  rsq = stats_trn$Rsq, rmse = stats_trn$RMSE,
#                  rmse.se = stats_trn$RMSE.se,
#                  bias = pbias(ta_prd, data_train$Ta_obs))
#     }))
# 
#     ndvimax <- stats[which.min(stats$rmse), "NDVImax"]
# 
#     Rsq_train <- stats[which.min(stats$rmse), "rsq"]
#     Rmse_train <- stats[which.min(stats$rmse), "rmse"]
#     Rmse_se_train <- stats[which.min(stats$rmse), "rmse.se"]
#     Bias_train <- stats[which.min(stats$rmse), "bias"]
# 
#     ## predict temperature
#     pred <- data_test$Slope * ndvimax + data_test$Intercept
#     resp <- data_test$Ta_obs
# 
#     stats_tst <- Rsenal::regressionStats(pred, resp, adj.rsq = FALSE)
# 
#     data.frame(n = n, NDVImax = ndvimax, TrainRMSE = Rmse_train,
#                TrainRMSEse = Rmse_se_train, TrainBias = Bias_train,
#                TrainRsq = Rsq_train, RMSE = stats_tst$RMSE,
#                RMSEse = stats_tst$RMSE.se, Bias = pbias(pred, resp),
#                Rsq = stats_tst$Rsq)
#   }))
# 
#   data.frame(habitat = h, t(data.frame(res)))
# }
# 
# saveRDS(dat_training, "data/results/stats_intra_7by7_area_terra.rds")
dat_training <- readRDS("data/results/stats_intra_7by7_area_terra.rds")

## reorder factor levels
dat_training$habitat <- factor(dat_training$habitat, 
                               levels = rev(sortElevation(FALSE)))


### evaluate approach at habitat scale -----

# ## select remaining plots
# dat_test <- dat_select[!dat_select$PlotID %in% dat_sel$PlotID, ]
# 
# prediction_stats <- foreach(i = as.character(dat_training$habitat),
#                             .combine = "rbind", .packages = lib) %dopar% {
# 
#   ndvimax <- dat_training$NDVImax[dat_training$habitat == i]
#   dat_test_sub <- subset(dat_test, habitat == i)
# 
#   dat_mrg <- foreach(h = as.character(dat_test_sub$PlotID),
#                      .combine = "rbind") %do% {
# 
#     ## extract daily maximum temperature
#     dat_plots %>%
#       filter(plotID == h & !is.na(Ta_200)) %>%
#       group_by(Date = as.Date(datetime)) %>%
#       filter(length(Ta_200) == 24) %>%
#       summarise(Ta_200 = max(Ta_200, na.rm = TRUE)) %>%
#       data.frame() -> dat_sub
# 
#     ## synchronise with modis time steps, i.e. aggregate 8-day mean maximum values
#     dat_agg <- foreach(j = unique(lubridate::year(dat_sub$Date)),
#                        .combine = "rbind") %do% {
# 
#       # merge available dates with continuous daily time series
#       dts0 <- seq(as.Date(paste0(j, "-01-01")), as.Date(paste0(j, "-12-31")), 1)
#       mrg <- merge(data.frame(Date = dts0), dat_sub, all.x = TRUE)
# 
#       # aggregate 8-day intervals
#       mrg$Cuts <- cut(mrg$Date, "8 days")
#       mrg <- mrg[complete.cases(mrg), ]
# 
#       mrg %>%
#         group_by(Cuts) %>%
#         summarise(Ta_200 = round(mean(Ta_200, na.rm = TRUE), 2)) %>%
#         data.frame()
#     }
# 
#     ## extract intercept and slope values
#     dts_lst <- MODIS::extractDate(names(rst_lst))$inputLayerDates
#     avl <- dts_lst %in% strftime(dat_agg$Cuts, "%Y%j")
# 
#     metrics <- foreach(k = list(mat_intercept, mat_slope),
#                        j = list(rst_intercept, rst_slope)) %do% {
# 
#       mat <- k[, avl]
#       if (class(mat) == "numeric")
#         mat <- matrix(mat, ncol = 1)
# 
#       rst <- j[[which(avl)]]
# 
#       cell_id <- cellFromXY(j, plots[plots@data$PlotID == h, ])
#       val <- mat[cell_id, ]
# 
#       if (all(is.na(val))) {
#         data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates,
#                    rep(NA, length(MODIS::extractDate(names(rst))$inputLayerDates)))
#       } else {
#         data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates, val)
#       }
#                        }
# 
#     metrics <- Reduce(function(...) merge(..., by = "Date", all = TRUE), metrics)
#     names(metrics)[2:3] <- c("Intercept", "Slope")
# 
#     ## merge and return
#     dat_agg$Date <- strftime(dat_agg$Cuts, "%Y%j")
#     merge(dat_agg, metrics, all.y = TRUE)
#   }
# 
#   n <- sum(complete.cases(dat_mrg))
# 
#   ta_obs <- dat_mrg$Ta_200
# 
#   # estimate temperature
#   ta_prd <- (dat_mrg$Slope * ndvimax + dat_mrg$Intercept)
# 
#   # calculate r-squared, rmse, and bias
#   stats <- Rsenal::regressionStats(ta_prd, ta_obs, adj.rsq = FALSE)
# 
#   data.frame(habitat = i, NDVImax = ndvimax,
#              rsq = stats$Rsq, rmse = stats$RMSE, rmse.se = stats$RMSE.se,
#              bias = pbias(ta_prd, ta_obs), n = n)
# }
# 
# saveRDS(prediction_stats, "data/results/stats_inter_7by7_area_terra.rds")
prediction_stats <- readRDS("data/results/stats_inter_7by7_terra.rds")

## reorder factor levels
prediction_stats$habitat <- factor(prediction_stats$habitat, 
                                   levels = rev(sortElevation(FALSE)))


### evaluate approach at plot scale -----

## select remaining plots
# prediction_stats_plot <- foreach(i = as.character(dat_training$habitat), 
#                                  .combine = "rbind", .packages = lib) %dopar% {
#                               
#   ndvimax <- dat_training$NDVImax[dat_training$habitat == i]
#   dat_test_sub <- subset(dat_test, habitat == i)
#   
#   foreach(h = as.character(dat_test_sub$PlotID), .combine = "rbind") %do% {
#                                                    
#     ## extract daily maximum temperature
#     dat_plots %>%
#       filter(plotID == h & !is.na(Ta_200)) %>%
#       group_by(Date = as.Date(datetime)) %>%
#       filter(length(Ta_200) == 24) %>%
#       summarise(Ta_200 = max(Ta_200, na.rm = TRUE)) %>%
#       data.frame() -> dat_sub
#                        
#     ## synchronise with modis time steps, i.e. aggregate 8-day mean maximum values 
#     dat_agg <- foreach(j = unique(lubridate::year(dat_sub$Date)), 
#                        .combine = "rbind") %do% {
#                          
#       # merge available dates with continuous daily time series
#       dts0 <- seq(as.Date(paste0(j, "-01-01")), as.Date(paste0(j, "-12-31")), 1)
#       mrg <- merge(data.frame(Date = dts0), dat_sub, all.x = TRUE)
#       
#       # aggregate 8-day intervals
#       mrg$Cuts <- cut(mrg$Date, "8 days")
#       mrg <- mrg[complete.cases(mrg), ]
#       
#       mrg %>% 
#         group_by(Cuts) %>%
#         summarise(Ta_200 = round(mean(Ta_200, na.rm = TRUE), 2)) %>%
#         data.frame()
#     }
#     
#     ## extract intercept and slope values
#     dts_lst <- MODIS::extractDate(names(rst_lst))$inputLayerDates
#     avl <- dts_lst %in% strftime(dat_agg$Cuts, "%Y%j")
#     
#     if (all(!avl) | sum(avl) == 1) {
#       return(data.frame(PlotID = h, NDVImax = NA, rsq = NA, rmse = NA, 
#                         stde = NA, bias = NA, n = ifelse(sum(avl) == 1, 1, NA)))
#     }
#     
#     metrics <- foreach(k = list(mat_intercept, mat_slope), 
#                        j = list(rst_intercept, rst_slope)) %do% {
#                            
#       mat <- k[, avl]
#       if (class(mat) == "numeric")
#         mat <- matrix(mat, ncol = 1)
#       
#       rst <- j[[which(avl)]]
#       
#       cell_id <- cellFromXY(j, plots[plots@data$PlotID == h, ])
#       val <- mat[cell_id, ]
#       
#       if (all(is.na(val))) {
#         data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates, 
#                    rep(NA, length(MODIS::extractDate(names(rst))$inputLayerDates)))
#       } else {
#         data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates, val)
#       }
#     }
#     
#     metrics <- Reduce(function(...) merge(..., by = "Date", all = TRUE), metrics)
#     names(metrics)[2:3] <- c("Intercept", "Slope")
#     
#     if (all(is.na(metrics$Intercept)) | all(is.na(metrics$Slope))) {
#       return(data.frame(PlotID = h, NDVImax = NA, rsq = NA, rmse = NA, 
#                         stde = NA, bias = NA, n = NA))
#     }
#     
#     ## merge and return
#     dat_agg$Date <- strftime(dat_agg$Cuts, "%Y%j")
#     dat_mrg <- merge(dat_agg, metrics, all.y = TRUE)
#     
#     n <- sum(complete.cases(dat_mrg))
#     
#     ta_obs <- dat_mrg$Ta_200
#     
#     # estimate temperature
#     ta_prd <- (dat_mrg$Slope * ndvimax + dat_mrg$Intercept)
#     
#     # calculate r-squared, rmse, and bias
#     data.frame(PlotID = h, NDVImax = ndvimax, 
#                rsq = summary(lm(ta_prd ~ ta_obs))$r.squared, 
#                rmse = rmse(ta_obs, ta_prd), stde = stde(ta_prd),
#                bias = pbias(ta_prd, ta_obs), n = n)
#   }
# }
#   
# # saveRDS(prediction_stats_plot, "data/results/stats_plot_inter_7by7_terra.rds")
# prediction_stats_plot <- readRDS("data/results/stats_plot_inter_7by7_terra.rds")
# 
# prediction_stats_plot <- merge(prediction_stats_plot, 
#                                dat_test[, c("PlotID", "habitat")], 
#                                by = "PlotID", all.x = TRUE)
# 
# ## reorder factor levels
# prediction_stats_plot$habitat <- factor(prediction_stats_plot$habitat, 
#                                    levels = c("mai", "sav", "cof", "hom", "gra", 
#                                               "flm", "fod", "foc", "fpo", "fpd", 
#                                               "fer", "fed", "hel"))
