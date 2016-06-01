library(raster)
library(doParallel)
library(rgeos)

## parallelization
cl <- makeCluster(detectCores() - 1)
registerDoParallel()

## polygons: ndvi
rst_tmp <- raster("data/MCD09Q1.A2011001.sur_refl.tif")
rst_tmp[is.na(rst_tmp[])] <- 0
spy_ndvi <- rasterToPolygons(rst_tmp)

## polygons: lst
rst_tmp <- raster("data/MOD11A2.A2011001.LST_Day_1km.tif")
rst_tmp[is.na(rst_tmp[])] <- 0
spy_lst <- rasterToPolygons(rst_tmp)

## calculate proportionate area contributions of modis cmg cell per gimms cell
lst_areas <- foreach(i = 1:length(spy_lst), 
                     .packages = c("raster", "rgeos")) %dopar% {
  cutset <- gIntersects(spy_ndvi, spy_lst[i, ], byid = TRUE)
  cutset <- as(cutset, "logical")
  spy_ndvi_sub <- spy_ndvi[cutset, ]
  
  spy_ndvi_areas <- sapply(1:length(spy_ndvi_sub), function(j) {
    gArea(intersect(spy_ndvi_sub[j, ], spy_lst[i, ])) / 
      gArea(spy_lst[j, ])
  })
  
  data.frame(cell = which(cutset), area = spy_ndvi_areas)
}

saveRDS(lst_areas, "data/weighted_area.rds")
