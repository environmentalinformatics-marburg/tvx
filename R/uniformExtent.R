uniformExtent <- function(template = "data/reference_grid.tif", f = 5000L, 
                          outProj = NULL) {
  
  ## required packages
  lib <- c("raster", "rgdal", "GSODTools")
  jnk <- sapply(lib, function(x) library(x, character.only = TRUE))
  
  ## import reference grid
  rst_ref <- raster(template)
  
  ## import coordinates of kilimanjaro intl airport
  kia <- subset(gsodstations, STATION.NAME == "KILIMANJARO INTL")
  coordinates(kia) <- ~ LON + LAT; proj4string(kia) <- "+init=epsg:4326"
  kia <- spTransform(kia, CRS = CRS(projection(rst_ref)))
  
  ## add expansion factor
  ext <- extent(rst_ref)
  ymin(ext) <- ymin(kia) - f
  ymax(ext) <- ymax(ext) + f
  xmin(ext) <- xmin(ext) - f
  xmax(ext) <- xmax(ext) + f
  
  if (!is.null(outProj)) {
    spy <- as(ext, "SpatialPolygons")
    proj4string(spy) <- projection(rst_ref)
    spy <- spTransform(spy, CRS = CRS(outProj))
    ext <- extent(spy)
  }
  
  return(ext)
}