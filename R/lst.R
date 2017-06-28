### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("doParallel", "MODIS")
Orcs::loadPkgs(lib)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

## modis options
MODISoptions("/media/fdetsch/XChange/MODIS_ARC", 
             "/media/fdetsch/XChange/MODIS_ARC/PROCESSED", 
             outProj = "+init=epsg:32637")


### data download and preprocessing -----

## reference extent (bale)
stn <- shapefile("../../data/bale/ema/ema_stations.shp")
ext <- as(extent(stn), "SpatialPolygons"); proj4string(ext) <- proj4string(stn)
bff <- rgeos::gBuffer(ext, width = 1e3)

## reference extent (kili)
rst_ref <- raster("data/reference_grid.tif")

## extract relevant sds
tfs <- runGdal("M*D11A2", extent = bff, job = "balelst",
               collection = getCollection("MOD11A2", forceCheck = TRUE),
               SDSstring = "111011100011")


### data processing -----

## loop over single product
lst <- foreach(product = tfs, n = 1:length(tfs)) %do% {
  
  ## status message
  cat("Commencing with", names(tfs)[n], "processing...\n")

    
  ### quality control ----------------------------------------------------------
  ### discard cloudy pixels based on companion quality information ('QC_Day', 
  ### 'QC_Night')
  
  dir_prd <- paste0("data/", names(tfs)[n])
  if (!dir.exists(dir_prd)) dir.create(dir_prd)
  
  dir_qc <- paste0(dir_prd, "/qc")
  if (!dir.exists(dir_qc)) dir.create(dir_qc)
  
  ## perform quality check for day and night separately
  lst <- unlist(lapply(c(1, 4), function(i) unlist(sapply(product, "[[", i))))
  qcs <- unlist(lapply(c(2, 5), function(i) unlist(sapply(product, "[[", i))))
  
  ## loop over layers
  fls_qc <- paste0(dir_qc, "/", basename(lst))
  
  if (all(file.exists(fls_qc))) {
    cat("All quality-controlled", names(tfs)[n], "files exist, skipping iteration...\n")
    stack(fls_qc)
  } else {
    cat("Commencing with", names(tfs)[n], "quality control...\n")
    stack(foreach(i = 1:length(lst), .packages = "raster", 
                  .export = ls(envir = globalenv())) %dopar% {
      
      if (file.exists(fls_qc[i])) {
        raster(fls_qc[k])
      } else {
        r1 <- raster(lst[i]); r1 <- r1 * 0.02 - 273.15
        r2 <- raster(qcs[i])
        
        overlay(r1, r2, fun = function(x, y) {
          id <- sapply(y[], function(l) {
            bin <- Rsenal::number2binary(l, 8, TRUE)
            mandatory_qa <- substr(bin, 7, 8)
            data_quality <- substr(bin, 5, 6)
            
            # pixel produced, good quality
            if (mandatory_qa == "00" | data_quality == "00") {
              return(TRUE)
              
              # pixel produced, unreliable or unquantifiable quality
            } else if (mandatory_qa == "01" | data_quality == "01") {
              emis_error <- substr(bin, 3, 4) == "00"
              lst_error <- substr(bin, 1, 2) == "00"
              
              return(all(emis_error, lst_error))
              
              # pixel not produced due to cloud effects or other reasons
            } else {
              return(FALSE)
            }
          })
          
          x[!id] <- NA
          return(x)
        }, filename = fls_qc[i])
      }
    })
  }
}
  

### combined product -----

## remove missing dates in terra-modis from aqua-modis stacks
rmv <- which(!MODIS::extractDate(names(lst[[2]][[1]]))$inputLayerDates %in% 
               MODIS::extractDate(names(lst[[1]][[1]]))$inputLayerDates)

lst[[2]][[1]] <- lst[[2]][[1]][[-rmv]]
lst[[2]][[2]] <- lst[[2]][[2]][[-rmv]]

lst_fill <- foreach(i = append(lst[[1]], lst[[2]]), 
                    j = append(lst[[2]], lst[[1]]), 
                    l = append(lst[[1]][c(2, 1)], lst[[2]][c(2, 1)]), 
                    m = append(lst[[2]][c(2, 1)], lst[[1]][c(2, 1)])) %do% {
                      
                      dir_fill <- paste0(dirname(dirname(attr(i[[1]], "file")@name)), "/gf")
                      if (!dir.exists(dir_fill)) dir.create(dir_fill)
                      
                      fls_fill <- paste0(dir_fill, "/", names(i), ".tif")
                      
                      if (all(file.exists(fls_fill))) {
                        stack(fls_fill)
                      } else {                      
                        
                        mat_resp <- raster::as.matrix(i)
                        mat_pred1 <- raster::as.matrix(j)
                        mat_pred2 <- raster::as.matrix(l)
                        mat_pred3 <- raster::as.matrix(m)
                        
                        mat_fill <- foreach(k = 1:nrow(mat_resp), .combine = "rbind") %do% {
                          dat <- data.frame(y = mat_resp[k, ], x1 = mat_pred1[k, ], 
                                            x2 = mat_pred2[k, ], x3 = mat_pred3[k, ])
                          
                          if (sum(complete.cases(dat)) >= (.25 * nrow(dat))) {
                            mod <- lm(y ~ x1 + x2 + x3, data = dat)
                            
                            id <- which(!is.na(dat$x1) & !is.na(dat$x2) & !is.na(dat$x3) & is.na(dat$y))
                            newdata <- data.frame(x1 = dat$x1[id], x2 = dat$x2[id], x3 = dat$x3[id])
                            dat$y[id] <- predict(mod, newdata)
                          }
                          
                          return(dat$y)
                        }
                        
                        rst_fill <- raster::setValues(i, mat_fill)
                        
                        lst_fill <- foreach(i = 1:ncol(mat_resp), .packages = "raster") %dopar%
                          raster::writeRaster(rst_fill[[i]], filename = fls_fill[i], 
                                              format = "GTiff", overwrite = TRUE)
                        
                        raster::stack(lst_fill)
                      }
                    }

## deregister parallel backend
stopCluster(cl)
