### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages and functions
lib <- c("doParallel", "gimms", "MODIS")
Orcs::loadPkgs(lib)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### download and rasterize data -----

## download data
fls <- downloadGimms(dsn = "../../programming/r/gimms/data", cores = 3L)

## reference extent from modis
ref <- raster("data/MCD13Q1.006/whittaker/MCD.2002177.yL5000.ndvi.tif")
ref <- projectExtent(ref, crs = projection(raster(fls[1], varname = "ndvi")))

## rasterize data and flags
drs_qc <- "data/GIMMS3g.v1/qc"
if (!dir.exists(drs_qc)) dir.create(drs_qc)

fls_ndvi <- paste(drs_qc, gsub("nc4$", "tif", basename(fls)), sep = "/")
rst_ndvi <- rasterizeGimms(fls, ext = ref, keep = 0, cores = 3L, 
                           filename = fls_ndvi, format = "GTiff")


### monthly aggregation -----

rst_mvc <- monthlyComposite(rst_ndvi, monthlyIndices(fls), cores = 3L)

dts_mvc <- monthlyIndices(fls, timestamp = TRUE)
dts_mvc <- dts_mvc[seq(1, length(dts_mvc), 2)]


### whittaker smoothing -----

## target folder and files
drs_wht <- "data/GIMMS3g.v1/wht"
if (!dir.exists(drs_wht)) dir.create(drs_wht)

fls_qc <- gsub("/rst/ndvi/", "/qc/", fls_ndvi)

## apply whittaker smoother
lst_wht <- whittaker.raster(rst_mvc, timeInfo = orgTime(dts_mvc), 
                            outDirPath = drs_wht, overwrite = TRUE, 
                            format = "GTiff")

## write whittaker-smoothed images to disc
rst_wht <- stack(lst_wht)

fls_wht <- basename(oldNaming(fls, ".tif"))
fls_wht <- fls_wht[seq(1, length(fls_wht), 2)]
fls_wht <- paste(drs_wht, fls_wht, sep = "/")

lst_wht <- foreach(i = 1:nlayers(rst_wht), .packages = "raster") %dopar% {
  rst <- rst_wht[[i]]
  rst[rst[] > 1] <- NA
  
  writeRaster(rst, filename = fls_wht[i], format = "GTiff", overwrite = TRUE)
}

rst_wht <- stack(lst_wht); rm(lst_wht)

## remove deprecated whittaker-related files
fls_old <- list.files(drs_wht, pattern = "yL5000.ndvi.tif", full.names = TRUE)
jnk <- file.remove(fls_old)

## deregister parallel backend
stopCluster(cl)
