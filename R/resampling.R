### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("Orcs", "doParallel", "MODIS", "Rsenal")
Orcs::loadPkgs(lib)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


### preprocessing -----

## lst
fls_lst <- list.files("data/MOD11A2.006/qc", pattern = "Day_1km.tif$", 
                      full.names = TRUE)
rst_lst <- stack(fls_lst)
dts_lst <- extractDate(fls_lst)$inputLayerDates
mat_lst <- as.matrix(rst_lst)

## 1-km ndvi (confirmed that Whittaker smoothed data yields better regression 
## metrics than 'qc2' data)
fls_ndvi <- list.files("data/MCD13A2.006/whittaker", pattern = ".tif$", 
                       full.names = TRUE)
rst_ndvi <- stack(fls_ndvi)
dts_ndvi <- extractDate(fls_ndvi)$inputLayerDates

id <- which(!(dts_ndvi %in% dts_lst))
rst_ndvi <- rst_ndvi[[-id]]; dts_ndvi <- dts_ndvi[-id]
mat_ndvi <- as.matrix(rst_ndvi)

## 250-m ndvi
fls_ndvi250 <- list.files("data/MCD13Q1.006/whittaker", pattern = ".tif$", 
                          full.names = TRUE)
rst_ndvi250 <- stack(fls_ndvi250)
dts_ndvi250 <- extractDate(fls_ndvi250)$inputLayerDates

id <- which(!(dts_ndvi250 %in% dts_lst))
rst_ndvi250 <- rst_ndvi250[[-id]]; dts_ndvi250 <- dts_ndvi250[-id]
mat_ndvi250 <- as.matrix(rst_ndvi250)

## dem
rst_dem <- raster("data/DEM/DEM_ARC1960_30m_Hemp.tif")
rst_dem <- resample(rst_dem, rst_ndvi)
val_dem <- rst_dem[]

rst_dem250 <- raster("data/DEM/DEM_ARC1960_30m_Hemp.tif")
rst_dem250 <- resample(rst_dem250, rst_ndvi250)
val_dem250 <- rst_dem250[]


### processing -----

## slope and intercept per scene
dat_cfs <- foreach(j = 1:ncol(mat_lst), .combine = "rbind") %do% {
  x <- mat_ndvi[, j]; y <- mat_lst[, j]
  if (!(all(is.na(x)) | all(is.na(y)))) {
    mod <- lm(y ~ x + val_dem)
    int <- coef(mod)[1]; slp1 <- coef(mod)[2]; slp2 <- coef(mod)[3]
    rsq <- summary(mod)$r.squared; p <- pvalue(mod)
  } else {
    int <- slp1 <- slp2 <- rsq <- p <- NA
  }

  data.frame(scene = j, date = extractDate(names(rst_lst[[j]]))$inputLayerDates, 
             intercept = int, slope1 = slp1, slope2 = slp2, rsq = rsq, p = p)
}

## identify 250-m cells within each 1-km cell
# dat_cls <- data.frame(matrix(nrow = ncell(rst_ndvi), ncol = 16))
# 
# for (i in 1:ncell(rst_ndvi)) {
#   if (i %% 100 == 0)
#     cat("Now processing cell no.", i, "...\n")
#   
#   rst <- rst_ndvi[[1]]
#   rst[i] <- 1; rst[][-i] <- NA
#   shp <- rasterToPolygons(rst)
#   cls <- unlist(cellFromPolygon(rst_ndvi250, shp))
#   dat_cls[i, ] <- cls
# }
# 
# saveRDS(dat_cls, "data/cells250m.rds")
dat_cls <- readRDS("data/cells250m.rds")

## perform resampling
dir_lst250 <- "data/MOD11Q1.006"
if (!dir.exists(dir_lst250)) dir.create(dir_lst250)

fls_lst250 <- paste0(dir_lst250, "/", basename(fls_lst))

lst_lst250 <- foreach(i = 1:nlayers(rst_lst), .packages = "raster") %dopar% {
  if (file.exists(fls_lst250[i])) {
    raster(fls_lst250[i])
  } else {
    rst <- rst_ndvi250[[i]] * dat_cfs$slope1[i] + rst_dem250 * dat_cfs$slope2[i] + 
      dat_cfs$intercept[i]
  
    writeRaster(rst, fls_lst250[i], format = "GTiff", overwrite = TRUE)
  }
}

rst_lst250 <- stack(lst_lst250); rm(lst_lst250)

## deregister parallel backend
stopCluster(cl)
