### environmental stuff -----

## load packages
Orcs::loadPkgs(c("raster", "rgdal", "doParallel", "dplyr"))

## load functions
source("R/movingModel.R")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### lst -----

## import files (daytime)
fls_lst <- list.files("data/MYD11A2.005/qc", 
                      pattern = "^MYD11A2.*_Day_1km.tif$", full.names = TRUE)

fls_lst <- fls_lst[grep("2011", fls_lst)[1]:grep("2015361", fls_lst)]

rst_lst <- stack(fls_lst) - 273.15
dts_lst <- MODIS::extractDate(fls_lst)$inputLayerDates

# ## import files (nighttime)
# fls_lst_night <- list.files("data/MYD11A2.005/gf", full.names = TRUE,
#                             pattern = "^MYD11A2.*_Night_1km.tif$")
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
rst_ndvi_res <- rst_ndvi_res / 10e3


# ### correlation -----
# 
# ## convert ndvi layers to matrix
# mat_ndvi_res <- raster::as.matrix(rst_ndvi_res)
# 
# ## loop over lst daytime and nighttime scenes
# lst_r <- foreach(i = list(rst_lst, rst_lst_night), 
#                  .packages = c("raster", "foreach", "Orcs")) %dopar% {
#                    
#   # convert current lst layers to matrix                 
#   mat_lst <- raster::as.matrix(i)
#   
#   # loop over matrix rows (i.e., pixels)
#   mat_r <- foreach(k = 1:nrow(mat_lst), .combine = "rbind") %do% {
#     
#     # if enough valid values are available, compute linear regression metrics
#     dat <- data.frame(y = mat_lst[k, ], x = mat_ndvi_res[k, ])
#     
#     if (sum(complete.cases(dat)) >= (.25 * nrow(dat))) {
#       dat <- dat[complete.cases(dat), ]
#       
#       r <- cor(dat$x, dat$y)
#       mod <- lm(y ~ x, data = dat)
#       
#       # return metrics only if highly significant (p < 0.001)
#       ifelse(Orcs::pvalue(mod) < .001, r, NA)
#       
#     } else {
#       return(NA)
#     }
#   }
#   
#   # write r values to raster
#   rst_r <- setValues(i[[1]], mat_r)
#   names(rst_r) <- ""
#   
#   return(rst_r)
# }
  

### temperature-vegetation index (tvx) related slope, intercept, and
### regression coefficient -----

## target folder
dir_tvx <- "data/MYD11A2.005/tvx"
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
    
    mod <- movingModel(rst_lst[[i]], rst_ndvi_res[[i]])
    mod <- stack(mod)
    
    # save metrics
    writeRaster(mod, filename = fls_tvx, format = "GTiff", overwrite = TRUE)
  }
}

          
### station data -----
          
## import plot coordinates
plots <- readOGR(dsn = "data/station_data", 
                 layer = "PlotPoles_ARC1960_mod_20140807_final", 
                 p4s = "+init=epsg:21037")

plots <- subset(plots, PoleType == "AMP")

## import available plot data 
dat <- read.csv("data/station_data/plots.csv")
dat$Ta_200 <- dat$Ta_200

# dat$datetime <- strptime(dat$datetime, format = "%Y-%m-%dT%H:%S")
# 
# dat_sub <- subset(dat, !is.na(Ta_200) & 
#   datetime == as.POSIXlt("2014-07-04 13:00", format = "%Y-%m-%d %H:%S"))
# 
# ## reorder plots with available dat 
# id <- do.call("c", lapply(dat_sub$plotID, function(i) {
#   ifelse(length(grep(i, plots$PlotID)) > 0, grep(i, plots$PlotID), NA)
# }))
# 
# if (any(is.na(id)))
#   dat_sub <- dat_sub[!is.na(id), ]
# 
# plots_sub <- plots[na.omit(id), ]


### evaluation -----

## sample station data: mai0
out <- lapply(as.character(unique(dat$plotID)), function(h) {
  cat("Processing plot", h, "...\n")
  
  dat %>%
    filter(plotID == h & !is.na(Ta_200)) %>%
    group_by(Date = as.Date(datetime)) %>%
    filter(length(Ta_200) == 24) %>%
    summarise(Ta_200 = max(Ta_200, na.rm = TRUE)) %>%
    data.frame() -> dat_sub
  
  ## synchronise with modis time steps, i.e. aggregate 8-day mean values 
  dat_agg <- foreach(i = unique(lubridate::year(dat_sub$Date)), 
                     .combine = "rbind") %do% {
                       
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
    return(data.frame(PlotID = h, NDVImax = NA, rsq = NA, rmse = NA))
  }
  
  metrics <- foreach(i = 2:1) %do% {
    rst <- stack(sapply(lst_tvx, "[[", i))
    rst <- rst[[which(avl)]]
    
    val <- as.numeric(extract(rst, plots[plots@data$PlotID == h, ]))
    
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
    return(data.frame(PlotID = h, NDVImax = NA, rsq = NA, rmse = NA))
  }
  
  ## estimate temperature using a priori ndvi_max (taken from czajkowski et al., 
  ## 2000; stisen et al., 2007)
  dat_agg$Date <- strftime(dat_agg$Cuts, "%Y%j")
  dat_mrg <- merge(dat_agg, metrics, all.y = TRUE)
  
  ta_obs <- dat_mrg$Ta_200
  
  stats <- lapply(seq(.2, 1, .01), function(i) {
    ta_prd <- (dat_mrg$Slope * i + dat_mrg$Intercept)
    
    df <- data.frame(PlotID = h, NDVImax = i, 
                      rsq = summary(lm(ta_prd ~ ta_obs))$r.squared, 
                      rmse = rmse(ta_obs, ta_prd))
    
    p <- xyplot(ta_prd ~ ta_obs, 
           xlab = expression("Ta"[obs]), ylab = expression("Ta"[prd]), 
           panel = function(x, y, ...) {
             panel.xyplot(x, y, col = "grey50", ...)
             panel.ablineq(lm(y ~ x), rotate = TRUE, r.squared = TRUE)
           })
    
    list(df, p)
  })
  
  id <- which.min(do.call("rbind", lapply(stats, "[[", 1))$rmse)
  stats[[id]]
})

saveRDS(out, file = "data/results/ndvimax_rmse.rds")

  
  ### test -----

dat_mrg %>% 
  mutate(NDVImax = (Ta_200 - Intercept) / Slope) %>%
  data.frame() -> dat_tmp


dat_tmp[dat_tmp$NDVImax > 1 & !is.na(dat_tmp$NDVImax), "NDVImax"] <- 1


rmse <- function(x, y) sqrt(mean((y - x)^2, na.rm = TRUE))

ndviMax(a, b, ta)

lsfit(x = rep(ndvimax0, length(ta)), y = ta)
lsfit(x = ta, y = rep(ndvimax0, length(ta)))

ndviMax(20.459, 0.005, ta)
