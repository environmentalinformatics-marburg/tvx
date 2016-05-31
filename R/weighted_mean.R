library(rgeos)

## polygons: ndvi
rst_tmp <- rst_ndvi[[1]]
rst_tmp[is.na(rst_tmp[])] <- 0
spy_ndvi <- rasterToPolygons(rst_tmp)

## polygons: lst
rst_tmp <- rst_lst[[1]]
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