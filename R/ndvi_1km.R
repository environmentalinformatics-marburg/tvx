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
MODISoptions(localArcPath = "/media/dogbert/modis_data/MODIS_ARC", 
             outDirPath = "/media/dogbert/modis_data/MODIS_ARC/PROCESSED", 
             outProj = "+init=epsg:21037")


### preprocessing -----

## download .hdf files in parallel
foreach(product = c("MOD13A2", "MYD13A2"), .packages = "MODIS") %dopar% {
  getHdf(product, tileH = 21, tileV = 9, begin = "2011001",
         collection = getCollection(product, forceCheck = TRUE))
}

## extract required sds
for (product in c("MOD13A2", "MYD13A2")) {
  runGdal(product, tileH = 21, tileV = 9, begin = "2011001",
          collection = getCollection(product, forceCheck = TRUE), 
          job = paste0(product, ".006"), SDSstring = "101000000001")
}


### processing -----

## reference grid
rst_ref <- raster("data/reference_grid.tif")

## loop over products
lst_prd <- lapply(c("MOD13A2", "MYD13A2"), function(product) {
  
  dir_prd <- paste0("data/", product, ".006")
  if (!dir.exists(dir_prd)) dir.create(dir_prd)
  
  ## crop images
  rst_crp <- foreach(i = c("NDVI", "pixel_reliability", "VI_Quality"), 
                     .packages = lib, .export = ls(envir = globalenv())) %dopar% {                                      
                       
    # list and import available files
    fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/", product, ".006"),
                      pattern = paste0(i, ".tif$"), full.names = TRUE)
    rst <- raster::stack(fls)
    
    # crop
    dir_out <- paste0(dir_prd, "/crp")
    if (!dir.exists(dir_out)) dir.create(dir_out)
    
    fls_out <- paste0(dir_out, "/", basename(fls))
    
    lst_out <- lapply(1:(raster::nlayers(rst)), function(j) {
      if (file.exists(fls_out[j])) {
        raster::raster(fls_out[j])
      } else {
        rst_out <- raster::crop(rst[[j]], rst_ref, snap = "out")
        
        # apply scale factor
        if (i == "NDVI")
          rst_out <- rst_out * 0.0001
        
        # save and return cropped layers
        raster::writeRaster(rst_out, filename = fls_out[j],
                            format = "GTiff", overwrite = TRUE)
      }
    })
    
    raster::stack(lst_out)
  }
  
  
  ### quality control, step #1: -----
  ### discard clouds, snow/ice and filled pixels using 'pixel_reliability'
  
  dir_qc1 <- paste0(dir_prd, "/qc1")
  if (!dir.exists(dir_qc1)) dir.create(dir_qc1)
  
  fls_qc1 <- paste0(dir_qc1, "/", names(rst_crp[[1]]), ".tif")
  
  ## perform quality check #1
  lst_qc1 <- foreach(i = 1:nlayers(rst_crp[[1]]), .packages = lib, 
                     .export = ls(envir = globalenv())) %dopar% {
    if (file.exists(fls_qc1[i])) {
      raster(fls_qc1[i])
    } else {
      overlay(rst_crp[[1]][[i]], rst_crp[[2]][[i]], fun = function(x, y) {
        x[!y[] %in% c(0, 1)] <- NA
        return(x)
      }, filename = fls_qc1[i], overwrite = TRUE, format = "GTiff")
    }
  }
  
  rst_qc1 <- stack(lst_qc1)
  
  
  ### quality control, step #2: -----
  ### discard cloudy pixels based on 'state_250m' flags
  
  dir_qc2 <- paste0(dir_prd, "/qc2")
  if (!dir.exists(dir_qc2)) dir.create(dir_qc2)
  
  fls_qc2 <- paste0(dir_qc2, "/", names(rst_qc1), ".tif")
  
  ## perform quality check #2
  lst_qc2 <- foreach(i = 1:nlayers(rst_qc1), .packages = lib, 
                     .export = ls(envir = globalenv())) %dopar% {
    if (file.exists(fls_qc2[i])) {
      raster(fls_qc2[i])
    } else {
      overlay(rst_qc1[[i]], rst_crp[[3]][[i]], fun = function(x, y) {
        id <- sapply(y[], function(k) {
          bin <- number2binary(k, 16, TRUE)
          quality <- substr(bin, 15, 16)
          
          if (quality == "00") {
            return(TRUE)
          } else if (quality %in% c("10", "11")) {
            return(FALSE)
          } else {
            useful <- !substr(bin, 11, 14) %in% c("1101", "1110")
            aerosol <- substr(bin, 9, 10) != "11"
            adjacent <- substr(bin, 8, 8) == "0"
            mixed <- substr(bin, 6, 6) == "0"
            snow <- substr(bin, 2, 2) == "0"
            shadow <- substr(bin, 1, 1) == "0"
            
            all(useful, aerosol, adjacent, mixed, snow, shadow)
          }
        })
        
        x[!id] <- NA
        return(x)
      }, filename = fls_qc2[i], overwrite = TRUE, format = "GTiff")
    }
  }
  
  raster::stack(lst_qc2)
})


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

## deregister parallel backend
stopCluster(cl)
