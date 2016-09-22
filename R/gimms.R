### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("gimms")
Orcs::loadPkgs(lib)


### processing -----

## download data
fls <- downloadGimms(dsn = "data/gimms3g/raw", cores = 3L)

## rasterize data and flags
fls_ndvi <- gsub("raw", "rst/ndvi", fls)
fls_ndvi <- paste0(fls_ndvi, ".tif")
rst_ndvi <- rasterizeGimms(fls, cores = 3L, filename = fls_ndvi)

fls_flag <- gsub("raw", "rst/flag", fls)
fls_flag <- paste0(fls_flag, ".tif")
rst_flag <- rasterizeGimms(fls, cores = 3L, filename = fls_flag)