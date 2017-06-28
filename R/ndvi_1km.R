### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("doParallel", "MODIS", "Rsenal")
Orcs::loadPkgs(lib)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## modis options
MODISoptions("/media/fdetsch/XChange/MODIS_ARC", 
             "/media/fdetsch/XChange/MODIS_ARC/PROCESSED", 
             outProj = "+init=epsg:32637")


### preprocessing -----

## reference extent (bale)
stn <- shapefile("../../data/bale/ema/ema_stations.shp")
ext <- as(extent(stn), "SpatialPolygons"); proj4string(ext) <- proj4string(stn)
bff <- rgeos::gBuffer(ext, width = 1e3)

## reference extent (kili)
rst_ref <- raster("data/reference_grid.tif")

## extract required sds
tfs <- runGdal("M*D13A2", extent = bff, job = "balendvi",
               collection = getCollection("MOD13A2", forceCheck = TRUE), 
               SDSstring = "101000000011")


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
