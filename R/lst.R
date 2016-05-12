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
MODISoptions(localArcPath = "data/MODIS_ARC", 
             outDirPath = "data/MODIS_ARC/PROCESSED", 
             outProj = "+init=epsg:21037")


### data download and preprocessing -----

# ## download data
# foreach(product = c("MOD11A2", "MYD11A2"), .packages = lib) %dopar%
#   MODIS::getHdf(product, tileH = 21, tileV = 9, collection = "005", 
#                 begin = "2011001")
# 
# ## extract relevant sds
# foreach(product = c("MOD11A2", "MYD11A2")) %do%
#   MODIS::runGdal(product, tileH = 21, tileV = 9, collection = "005", 
#                  begin = "2011001",
#                  SDSstring = "110011000011", job = paste0(product, ".005"))


### data processing -----

## reference grid
rst_ref <- raster("data/reference_grid.tif")

## loop over single product
lst <- lapply(c("MOD11A2", "MYD11A2"), function(product) {
  
  ## status message
  cat("Commencing with the processing of", product, "...\n")
  
  ### crop layers ----------------------------------------------------------------
  
  ## setup output folder
  dir_out <- paste0("data/", product, ".005")
  if (!dir.exists(dir_out)) dir.create(dir_out)
  
  ## perform crop
  pattern <- c("Day_1km", "QC_Day", "Night_1km", "QC_Night")
  
  # rst_crp <- foreach(i = pattern, .packages = "MODIS") %dopar% {
  # 
  #   # list and import available files
  #   fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/", product, ".005"),
  #                     pattern = paste0(i, ".tif$"), full.names = TRUE)
  #   
  #   rst <- raster::stack(fls)
  # 
  #   # crop
  #   dir_crp <- paste0(dir_out, "/crp")
  #   if (!dir.exists(dir_crp)) dir.create(dir_crp)
  # 
  #   fls_crp <- paste(dir_crp, basename(fls), sep = "/")
  #   rst_crp <- raster::crop(rst, rst_ref, snap = "out")
  # 
  #   # if dealing with (day or night) lst bands, convert to 16-bit unsigned
  #   # integer and apply scale factor of 0.02
  #   if (i %in% c("Day_1km", "Night_1km")) {
  #     raster::dataType(rst_crp) <- "INT2U"
  #     rst_crp <- rst_crp * 0.02
  # 
  #   # else convert to 8-bit unsigned integer
  #   } else {
  #     raster::dataType(rst_crp) <- "INT1U"
  #   }
  # 
  #   # save and return cropped layers
  #   lst_crp <- lapply(1:nlayers(rst_crp), function(j)
  #     raster::writeRaster(rst_crp[[j]], filename = fls_crp[j],
  #                         format = "GTiff", overwrite = TRUE)
  #   )
  # 
  #   raster::stack(lst_crp)
  # }
  
  # ## reimport cropped files
  # dir_crp <- paste0(dir_out, "/crp")
  # 
  # rst_crp <- foreach(i = pattern, .packages = "raster") %dopar% {
  #   fls_crp <- list.files(dir_crp, pattern = paste0(i, ".tif$"), full.names = TRUE)
  #   raster::stack(fls_crp)
  # }
  
  
  ### quality control ----------------------------------------------------------
  ### discard cloudy pixels based on companion quality information ('QC_Day', 
  ### 'QC_Night')
  
  dir_qc <- paste0(dir_out, "/qc")
  if (!dir.exists(dir_qc)) dir.create(dir_qc)
  
  # ## perform quality check for day and night separately
  # lst_qc <- foreach(i = rst_crp[c(1, 3)], j = rst_crp[c(2, 4)]) %do% {
  # 
  #   ## loop over layers
  #   lst_out <- foreach(k = 1:nlayers(i), .packages = "raster") %dopar%
  #     overlay(i[[k]], j[[k]], fun = function(x, y) {
  #       id <- sapply(y[], function(l) {
  #         bin <- Rsenal::number2binary(l, 8, TRUE)
  #         mandatory_qa <- substr(bin, 7, 8)
  #         data_quality <- substr(bin, 5, 6)
  # 
  #         # pixel produced, good quality
  #         if (mandatory_qa == "00" | data_quality == "00") {
  #           return(TRUE)
  # 
  #         # pixel produced, unreliable or unquantifiable quality
  #         } else if (mandatory_qa == "01" | data_quality == "01") {
  #           emis_error <- substr(bin, 3, 4) == "00"
  #           lst_error <- substr(bin, 1, 2) == "00"
  # 
  #           return(all(emis_error, lst_error))
  # 
  #         # pixel not produced due to cloud effects or other reasons
  #         } else {
  #           return(FALSE)
  #         }
  #       })
  # 
  #       x[!id] <- NA
  #       return(x)
  #     }, filename = paste(dir_qc, names(i[[k]]), sep = "/"),
  #     overwrite = TRUE, format = "GTiff")
  # 
  #   raster::stack(lst_out)
  # }
  
  ## reimport quality-controlled files
  lst_qc <- lapply(pattern[c(1, 3)], function(i) {
    fls_qc <- list.files(dir_qc, pattern = paste0(i, ".tif$"), full.names = TRUE)
    fls_qc <- fls_qc[grep("2002185", fls_qc):grep("2016105", fls_qc)]
    raster::stack(fls_qc)
  })
  
  return(lst_qc)
})


### combined product -----------------------------------------------------------

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
