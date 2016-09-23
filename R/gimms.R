### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages and functions
lib <- c("doParallel", "gimms", "MODIS")
Orcs::loadPkgs(lib)

source("R/uniformExtent.R")

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### download and rasterize data -----

## download data
fls <- downloadGimms(dsn = "data/gimms3g/raw", cores = 3L)

## rasterize data and flags
fls_ndvi <- gsub("raw", "rst/ndvi", fls)
fls_ndvi <- paste0(fls_ndvi, ".tif")
rst_ndvi <- if (all(file.exists(fls_ndvi))) {
  stack(fls_ndvi)
} else {
  rasterizeGimms(fls, cores = 3L, filename = fls_ndvi)
}

fls_flag <- gsub("raw", "rst/flag", fls)
fls_flag <- paste0(fls_flag, ".tif")
rst_flag <- if (all(file.exists(fls_flag))) {
  stack(fls_flag)
} else {
  rasterizeGimms(fls, cores = 3L, filename = fls_flag)
}


### crop -----

## reference extent
ext <- uniformExtent(outProj = "+init=epsg:4326")

## loop over ndvi and flags
dir_crp <- "data/gimms3g/crp"
if (!dir.exists(dir_crp)) dir.create(dir_crp)

lst_crp <- foreach(i = list(rst_ndvi, rst_flag), j = list(fls_ndvi, fls_flag), 
                   .packages = lib) %dopar% {
                       
  # target folder and files
  fls_crp <- gsub("/rst/", "/crp/", j)
  if (!dir.exists(unique(dirname(fls_crp))))
    dir.create(unique(dirname(fls_crp)))
  
  # crop images
  lst_crp <- lapply(1:(raster::nlayers(i)), function(k) {
    if (file.exists(fls_crp[k])) {
      raster::raster(fls_crp[k])
    } else {
      raster::crop(i[[k]], ext, snap = "out", filename = fls_crp[k],
                   format = "GTiff", overwrite = TRUE)
    }
  })
  
  raster::stack(lst_crp)
}


### quality control -----

## target folder and files
dir_qc <- "data/gimms3g/qc"
if (!dir.exists(dir_qc)) dir.create(dir_qc)

fls_qc <- gsub("/rst/ndvi/", "/qc/", fls_ndvi)

## overlay ndvi with flags and keep 'good' values only (flags 1 and 2)
lst_qc <- foreach(i = 1:nlayers(lst_crp[[1]]), .packages = lib) %dopar% {
  if (file.exists(fls_qc[i])) {
    raster(fls_qc[i])
  } else {
    overlay(lst_crp[[1]][[i]], lst_crp[[2]][[i]], fun = function(x, y) {
      x[y[] > 2] <- NA
      return(x)
    }, filename = fls_qc[i], format = "GTiff", overwrite = TRUE)
  }
}

rst_qc <- stack(lst_qc); rm(lst_qc)


### whittaker smoothing --------------------------------------------------------

## set system locale
systime_locale <- Sys.getlocale(category = "LC_TIME")
if (Sys.info()[["sysname"]] == "Windows") {
  invisible(Sys.setlocale(category = "LC_TIME", locale = "C"))
} else {
  invisible(Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8"))
}

## replace %Y%m with %Y%m%d (compatible to `as.Date` in `orgTime`)
org_gimms <- basename(fls_ndvi)
for (i in 1:length(org_gimms)) {
  dt_yrmn <- substr(org_gimms[i], 4, 8)
  dt_yrmndy <- paste0(dt_yrmn, ifelse(i %% 2 == 1, "01", "15"))
  org_gimms[i] <- gsub(dt_yrmn, dt_yrmndy, org_gimms[i])
}

org_gimms <- orgTime(org_gimms, pillow = 0, pos1 = 4, pos2 = 10, format = "%y%b%d")

## target folder and files
dir_wht <- "data/gimms3g/whittaker"
if (!dir.exists(dir_wht)) dir.create(dir_wht)

fls_qc <- gsub("/rst/ndvi/", "/qc/", fls_ndvi)

## apply whittaker smoother
lst_wht <- whittaker.raster(rst_qc, timeInfo = org_gimms, outDirPath = dir_wht,
                            overwrite = TRUE, format = "GTiff")

## write whittaker-smoothed images to disc
rst_wht <- stack(lst_wht)
names(rst_wht) <- names(rst_qc)
fls_wht <- paste0(dir_wht, "/", names(rst_wht), ".tif")

lst_wht <- foreach(i = 1:nlayers(rst_wht), .packages = "raster") %dopar% {
  rst <- rst_wht[[i]]
  rst[rst[] > 1] <- NA
  
  writeRaster(rst, filename = fls_wht[i], format = "GTiff", overwrite = TRUE)
}

rst_wht <- stack(lst_wht); rm(lst_wht)

## remove deprecated whittaker-related files
fls_old <- list.files(dir_wht, pattern = "yL5000.ndvi.tif", full.names = TRUE)
file.remove(fls_old)


### monthly aggregation --------------------------------------------------------

## target folder and files
dir_mvc <- paste0(dir_wht, "/mvc")
if (!dir.exists(dir_mvc)) dir.create(dir_mvc)

fls_mvc <- gsub("/whittaker/", "/whittaker/mvc/", fls_wht)
fls_mvc <- gsub("15a", "15", fls_mvc)
fls_mvc <- gsub("15b", "15", fls_mvc)

lst_mvc <- strsplit(fls_mvc, "\\.")
fls_mvc <- sapply(1:length(lst_mvc), function(i) {
  paste(lst_mvc[[i]][-2], collapse = ".")
})

fls_mvc <- unique(fls_mvc)

## create monthly maximum value composites
rst_mvc <- monthlyComposite(fls_wht, cores = 3L, filename = fls_mvc, 
                            format = "GTiff", overwrite = TRUE)

## deregister parallel backend
stopCluster(cl)
