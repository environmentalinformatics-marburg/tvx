### environment ----------------------------------------------------------------

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("MODIS")
Orcs::loadPkgs(lib)

## modis options
MODISoptions(localArcPath = "data/MODIS_ARC", 
             outDirPath = "data/MODIS_ARC/PROCESSED", 
             outProj = "+init=epsg:21037")


### data download and preprocessing -----

collection <- getCollection("MCD12Q1", forceCheck = TRUE)

runGdal("MCD12Q1", collection, begin = "2013001", tileH = 21, tileV = 9, 
        job = "MCD12Q1.051", 
        SDSstring = paste(c(rep(1, 5), rep(0, 11)), collapse = ""))
        
## reference grid
rst_ref <- raster("data/reference_grid.tif")

## setup output folder
dir_out <- "data/MCD12Q1.051"
if (!dir.exists(dir_out)) dir.create(dir_out)

## list and import available files
fls <- list.files(paste0(getOption("MODIS_outDirPath"), "/MCD12Q1.051"),
                  pattern = ".tif$", full.names = TRUE)

rst <- stack(fls)

## crop
dir_crp <- paste0(dir_out, "/crp")
if (!dir.exists(dir_crp)) dir.create(dir_crp)

fls_crp <- paste(dir_crp, basename(fls), sep = "/")
rst_crp <- crop(rst, rst_ref, snap = "out")

## chala
spy_chala <- readOGR("data/station_data", "Lake_Chala")
spy_chala <- spTransform(spy_chala, CRS = CRS("+init=epsg:21037"))

id <- unlist(cellFromPolygon(rst_crp, spy_chala))
rst_crp[id] <- 0

lapply(1:nlayers(rst_crp), function(i) {
  rst_crp <- writeRaster(rst_crp[[i]], 
                         paste0("data/MCD12Q1.051/", names(rst_crp)[i]), 
                         format = "GTiff", overwrite = TRUE)
})
