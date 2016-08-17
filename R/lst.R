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
MODISoptions(localArcPath = "/media/dogbert/modis_data/MODIS_ARC", 
             outDirPath = "/media/dogbert/modis_data/MODIS_ARC/PROCESSED", 
             outProj = "+init=epsg:21037")


### data download and preprocessing -----

## download data
foreach(product = c("MOD11A2", "MYD11A2"), .packages = lib) %dopar%
  MODIS::getHdf(product, tileH = 21, tileV = 9, begin = "2011001", 
                collection = getCollection(product, forceCheck = TRUE))

## extract relevant sds
foreach(product = c("MOD11A2", "MYD11A2")) %do%
  MODIS::runGdal(product, tileH = 21, tileV = 9, begin = "2011001", 
                 collection = getCollection(product, forceCheck = TRUE),
                 SDSstring = "110011000011", job = paste0(product, ".006"))


### data processing -----

## reference grid
rst_ref <- raster("data/reference_grid.tif")

## loop over single product
lst <- lapply(c("MOD11A2", "MYD11A2"), function(product) {
  
  ## status message
  cat("Commencing with the processing of", product, "...\n")
  
  ### crop layers ----------------------------------------------------------------
  
  ## setup output folder
  dir_out <- paste0("data/", product, ".006")
  if (!dir.exists(dir_out)) dir.create(dir_out)
  
  ## perform crop
  pattern <- c("Day_1km", "QC_Day", "Clear_sky_days", 
               "Night_1km", "QC_Night", "Clear_sky_nights")
  
  rst_crp <- foreach(i = pattern, .packages = "MODIS", 
                     .export = ls(envir = globalenv())) %dopar% {
    
    # list and import available files
    fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/", product, ".006"),
                      pattern = paste0(i, ".tif$"), full.names = TRUE)
    
    rst <- raster::stack(fls)
    
    # crop
    dir_crp <- paste0(dir_out, "/crp")
    if (!dir.exists(dir_crp)) dir.create(dir_crp)
    
    fls_crp <- paste(dir_crp, basename(fls), sep = "/")

    lst_crp <- lapply(1:(raster::nlayers(rst)), function(j) {
      if (file.exists(fls_crp[j])) {
        raster::raster(fls_crp[j])
      } else {
        rst_crp <- raster::crop(rst[[j]], rst_ref, snap = "near")
        
        # if dealing with (day or night) lst bands, convert to 16-bit unsigned
        # integer and apply scale factor of 0.02 and offset of -273.15
        if (i %in% c("Day_1km", "Night_1km")) {
          raster::dataType(rst_crp) <- "INT2U"
          rst_crp <- rst_crp * 0.02 - 273.15
        
        # else if dealing with (day or night) no. of clear-sky observations, 
        # convert to 16-bit unsigned integer and apply scale factor of 0.0005
        } else if (i %in% c("Clear_sky_days", "Clear_sky_nights")) {
          raster::dataType(rst_crp) <- "INT2U"
          rst_crp <- rst_crp * 0.0005

        # else convert to 8-bit unsigned integer
        } else {
          raster::dataType(rst_crp) <- "INT1U"
        }
        
        # save and return cropped layers
        raster::writeRaster(rst_crp, filename = fls_crp[j],
                            format = "GTiff", overwrite = TRUE)
      }
    })
    
    raster::stack(lst_crp)
  }
  

  ### quality control ----------------------------------------------------------
  ### discard cloudy pixels based on companion quality information ('QC_Day', 
  ### 'QC_Night')
  
  dir_qc <- paste0(dir_out, "/qc")
  if (!dir.exists(dir_qc)) dir.create(dir_qc)
  
  ## perform quality check for day and night separately
  lst_qc <- foreach(i = rst_crp[c(1, 4)], j = rst_crp[c(2, 6)]) %do% {

    ## loop over layers
    fls_qc <- paste0(dir_qc, "/", names(i), ".tif")
    lst_out <- foreach(k = 1:nlayers(i), .packages = "raster", 
                       .export = ls(envir = globalenv())) %dopar% {
      
      if (file.exists(fls_qc[k])) {
        raster(fls_qc[k])
      } else {
        overlay(i[[k]], j[[k]], fun = function(x, y) {
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
        }, filename = fls_qc[k], overwrite = TRUE, format = "GTiff")
      }
    }

    stack(lst_out)
  }
})


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
