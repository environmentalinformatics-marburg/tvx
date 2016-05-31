#' Setup Moving Window
#' 
movingWindow <- function(width = 7L, height = width) {
  mat <- matrix(1, nrow = height, ncol = width)
  mat[ceiling(height/2), ceiling(width/2)] <- 0
  mat
}

#' Sort Habitats by Elevation
#' 
sortElevation <- function(df = TRUE) {
  
  ## load required packages
  library(dplyr)
  
  ## import plot data
  shp_plots <- readOGR("data/station_data", p4s = "+init=epsg:21037",
                       layer = "PlotPoles_ARC1960_mod_20140807_final")
  
  shp_plots <- subset(shp_plots, PoleType == "AMP")
  
  ## subset analyzed plots
  dat_select <- readRDS("data/results/ndvimax_rmse_7by7_terra.rds")
  shp_plots <- shp_plots[shp_plots$PlotID %in% dat_select$PlotID, ]
  
  ## set land-cover types for mch0 and mwh0
  shp_plots$Habitat <- substr(shp_plots$PlotID, 1, 3)
  shp_plots$Habitat[shp_plots$PlotID == "mch0"] <- "fpo"
  shp_plots$Habitat[shp_plots$PlotID == "mwh0"] <- "fpd"
  
  ## retrieve altitudinal range
  shp_plots@data %>%
    group_by(Habitat) %>%
    summarise(zmin = round(min(Z_DEM_HMP, na.rm = TRUE)), 
              zmax = round(max(Z_DEM_HMP, na.rm = TRUE))) %>%
    arrange(desc(zmin)) %>% 
    mutate(Elevation = paste(zmin, zmax, sep = "--")) %>%
    dplyr::select(Habitat, Elevation) %>%
    data.frame() -> dat_ele
  
  ## return sorted data.frame or vector of sorted habitats
  if (df) return(dat_ele) else return(dat_ele$Habitat)
}