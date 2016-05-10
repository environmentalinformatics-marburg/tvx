## create panel.ablineq for all available plots
dts_lst <- MODIS::extractDate(names(rst_lst))$inputLayerDates

lst_p <- lapply(as.character(unique(dat$plotID)), function(h) {
  
  ## status message
  cat("Commencing with plot", h, "...\n")
  
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
    
    if (all(is.na(mrg$Ta_200))) return(NULL)

    # aggregate 8-day intervals
    mrg$Cuts <- cut(mrg$Date, "8 days")
    mrg <- mrg[complete.cases(mrg), ]
    
    mrg %>% 
      group_by(Cuts) %>%
      summarise(Ta_200 = round(mean(Ta_200, na.rm = TRUE), 2)) %>%
      data.frame()
  }
  
  ## extract intercept and slope values
  avl <- dts_lst %in% strftime(dat_agg$Cuts, "%Y%j")
  
  if (all(!avl)) return(NULL)
  
  metrics <- foreach(i = 2:1) %do% {
    rst <- stack(sapply(lst_tvx, "[[", i))
    rst <- rst[[which(avl)]]
    
    data.frame(Date = MODIS::extractDate(names(rst))$inputLayerDates, 
               as.numeric(extract(rst, plots[plots@data$PlotID == "sav5", ])))
  }
  
  metrics <- Reduce(function(...) merge(..., by = "Date", all = TRUE), metrics)
  names(metrics)[2:3] <- c("Intercept", "Slope")
  
  ## estimate temperature using a priori ndvi_max (taken from czajkowski et al., 
  ## 2000; stisen et al., 2007)
  dat_agg$Date <- strftime(dat_agg$Cuts, "%Y%j")
  dat_mrg <- merge(dat_agg, metrics, all.y = TRUE)
  
  ndvimax0 <- .65
  
  ta_prd <- (dat_mrg$Slope * ndvimax0 + dat_mrg$Intercept) - 273.15
  ta_obs <- dat_mrg$Ta_200
  
  xyplot(ta_prd ~ ta_obs, 
         xlab = expression("Ta"[obs]), ylab = expression("Ta"[prd]), 
         panel = function(x, y, ...) {
           panel.xyplot(x, y, col = "grey50", ...)
           panel.ablineq(lm(y ~ x), rotate = TRUE, r.squared = TRUE)
         })
})

## remove unavailable plots
dat <- dat[dat$plotID %in% unique(dat$plotID)[sapply(lst_p, class) == "trellis"], ]
lst_p <- lst_p[sapply(lst_p, class) == "trellis"]


## merge available plots with coordinates data
plots_mrg <- merge(data.frame(PlotID = unique(dat$plotID)), 
                   plots@data, by = "PlotID", all.x = TRUE)

plots_mrg <- plots_mrg[complete.cases(plots_mrg), ]

coordinates(plots_mrg) <- ~ X + Y
proj4string(plots_mrg) <- "+init=epsg:21037"

## display data
mapview(plots_mrg, popup = popupGraph(lst_p, type = "svg", width = 4, height = 3))


###

out <- readRDS("data/results/ndvimax_rmse_(terra).rds")

lst_format <- lapply(out, "[[", 1)
lst_format <- lst_format[sapply(lst_format, class) == "data.frame"]
dat_format <- do.call("rbind", lst_format)

lst_p <- lapply(out, "[[", 2)
lst_p <- lst_p[sapply(lst_p, class) == "trellis"]

dat <- data.frame(plots)
dat <- merge(dat_format, dat, all.x = TRUE, sort = FALSE)

coordinates(dat) <- ~ X + Y
proj4string(dat) <- "+init=epsg:21037"

## display data
mapview(dat, popup = popupGraph(lst_p, type = "svg", width = 4, height = 3))
