### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("doParallel", "MODIS", "Rsenal", "ESD")
Orcs::loadPkgs(lib)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## modis options
lap <- "/media/fdetsch/XChange/MODIS_ARC"
MODISoptions(lap, paste0(lap, "/PROCESSED"), quiet = TRUE)


### preprocessing -----

## reference extent (bale)
shp <- shapefile("../../data/bale/ema/ema_stations.shp")
ext <- as(extent(shp), "SpatialPolygons"); proj4string(ext) <- proj4string(shp)
ema <- rgeos::gBuffer(ext, width = 1e4)

## reference extent (kili)
rst <- raster("data/reference_grid.tif"); ext <- extent(rst)
ymin(ext) <- ymin(ext) - 1e4
tma <- as(ext, "SpatialPolygons"); proj4string(tma) <- projection(rst)

## extract required sds
tfs <- foreach(i = list(ema, tma), j = as.list(1:2)) %do%
  runGdal("M*D13A2", extent = i, collection = "006", SDSstring = "101000000011",
          job = paste0("ndvi-1km-", ifelse(j == 1, "bale", "kili")))


### processing -----

## scale, perform quality control, create monthly maximum value composites and 
## combined terra/aqua images
lst_prd <- preprocessMODIS(tfs, dsn = "data", interval = "month", whit = FALSE, 
                           cores = 3L)


### whittaker smoother -----

## target folders and files
dir_prd <- "data/MCD13A2.006"
if (!dir.exists(dir_prd)) dir.create(dir_prd)

dir_wht <- paste0(dir_prd, "/whittaker")
if (!dir.exists(dir_wht)) dir.create(dir_wht)

# ## select temporal range
# for (i in seq(lst_prd)) {
#   st <- grep("2011", names(lst_prd[[i]]))[1]
#   nd <- grep("2015", names(lst_prd[[i]])); nd <- nd[length(nd)]
#   lst_prd[[i]] <- lst_prd[[i]][[st:nd]]
# }

## reorder layers
nms_qc2 <- do.call("c", lapply(lst_prd, names))
dts_qc2 <- extractDate(nms_qc2)$inputLayerDates
rst_qc2 <- stack(lst_prd)
rst_qc2 <- rst_qc2[[order(dts_qc2)]]
nms_qc2 <- nms_qc2[order(dts_qc2)]

## start with start date of aqua availability
rst_qc2 <- rst_qc2[[(grep("MYD13A2", nms_qc2)[1]-1):nlayers(rst_qc2)]]
nms_qc2 <- nms_qc2[(grep("MYD13A2", nms_qc2)[1]-1):length(nms_qc2)]

detach("package:MODIS", unload = TRUE)
install.packages("inst/extdata/MODIS_0.10-18.tar.gz",
                 repos = NULL, type = "source")
library(MODIS)

## apply whittaker smoother
lst_wht <- whittaker.raster(rst_qc2, outDirPath = dir_wht,
                            overwrite = TRUE, format = "GTiff")

## write to disc
rst_wht <- stack(lst_wht)
names(rst_wht) <- gsub("MOD13A2", "MCD13A2", nms_qc2)
names(rst_wht) <- gsub("MYD13A2", "MCD13A2", names(rst_wht))
fls_wht <- paste0(dir_wht, "/", names(rst_wht), ".tif")

lst_wht <- foreach(i = 1:nlayers(rst_wht), .packages = "raster") %dopar% {
  rst <- rst_wht[[i]]
  rst[rst[] > 1] <- NA
  
  writeRaster(rst, filename = fls_wht[i], format = "GTiff", overwrite = TRUE)
}

rst_wht <- stack(lst_wht)

## remove deprecated whittaker-related files
fls_old <- list.files(dir_wht, pattern = "NDVI_YearlyLambda", 
                      full.names = TRUE)
file.remove(fls_old)

## re-install new modis version
detach("package:MODIS", unload = TRUE)
devtools::install_github("MatMatt/MODIS", ref = "develop")

## deregister parallel backend
stopCluster(cl)
