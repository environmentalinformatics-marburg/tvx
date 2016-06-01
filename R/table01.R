### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("rgdal", "dplyr", "foreach", "stargazer")
Orcs::loadPkgs(lib)

## load functions
source("R/funs.R")


### ndvi_max -----

dat_all <- readRDS("data/results/stats_intra_7by7_area_terra.rds")[, 1:3]
dat_all$NDVImax <- round(dat_all$NDVImax, 2)
# dat_all$"Training plots" <- c("cof1,3", "fed3,4", "fer0", "flm1,4", "foc1,2", 
#                               "fod3,5", "fpd3,5", "fpo4,5", "gra4,5", 
#                               "hel4,5", "hom4,5", "mai3,4", "sav1,3")
# 
# dat_all <- dat_all[, c(1, 3, 2)]
# rownames(dat_all) <- NULL
names(dat_all)[1] <- "Habitat"


### n (remote)

dat_rem <- readRDS("data/results/stats_inter_7by7_area_terra.rds")[, c(1, 7)]

dat_all <- merge(dat_all, dat_rem, by.x = "Habitat", by.y = "habitat", 
                 all = TRUE, sort = FALSE)
names(dat_all)[c(2, 4)] <- c("n", "n_eval")

### elevations -----

dat_ele <- sortElevation()
dat_all <- merge(dat_ele, dat_all, by = "Habitat", all = TRUE, sort = FALSE)


### disturbance status -----

## disturbed habitats
dst <- c("fed", "fpd", "fod", "gra", "hom", "cof", "mai")

dat_all$Disturbance <- "No"
dat_all$Disturbance[dat_all$Habitat %in% dst] <- "Yes"


### ndvi-lst relationship -----

## selected plots
dat_select <- readRDS("data/results/ndvimax_rmse_7by7_area_terra.rds")

dat_select$habitat <- substr(dat_select$PlotID, 1, 3)
dat_select$habitat[dat_select$PlotID == "mch0"] <- "fpo"
dat_select$habitat[dat_select$PlotID == "mwh0"] <- "fpd"

dat_select %>%
  filter(n >= 10) %>%
  mutate(ratio = rmse / n * (-1)) %>%
  group_by(habitat) %>%
  top_n(n = 2, wt = ratio) %>%
  data.frame() -> dat_sel

## merge with spatial information
spt_plots <- readOGR("data/station_data", 
                     "PlotPoles_ARC1960_mod_20140807_final", 
                     p4s = "+init=epsg:21037")

spt_plots <- subset(spt_plots, PoleType == "AMP")

spt_r <- merge(dat_sel, spt_plots, by = "PlotID")
coordinates(spt_r) <- ~ X + Y
proj4string(spt_r) <- "+init=epsg:21037"

## extract r values
rst_r <- raster("data/results/ndvi_lst_r.tif")
rst_p <- raster("data/results/ndvi_lst_p.tif")

spt_r$r <- extract(rst_r, spt_r)
spt_r@data %>%
  group_by(habitat) %>%
  summarise(r = round(mean(r), 2)) %>%
  data.frame() -> dat_r

spt_r$p <- extract(rst_p, spt_r)
spt_r@data %>%
  group_by(habitat) %>%
  summarise(p = round(mean(p), 4)) %>%
  data.frame() -> dat_p

dat_all <- merge(dat_all, dat_r, by.x = "Habitat", by.y = "habitat", 
                 sort = FALSE)

dat_all$r <- as.character(dat_all$r)
for (i in 1:nrow(dat_p)) {
  id <- grep(dat_p$habitat[i], dat_all$Habitat)
  
  if (dat_p$p[i] < 0.001) {
    dat_all$r[id] <- paste0(dat_all$r[id], "*")
  }
}


### long habitat names -----

hbt <- c("Helichrysum", "Erica forest disturbed", "Erica forest", 
         "Podocarpus forest disturbed", "Podocarpus forest", 
         "Ocotea forest disturbed", "Ocotea forest", 
         "Lower montane forest", "Grassland", "Chagga homegarden", 
         "Coffee plantation", "Savanna", "Maize field")

jnk <- foreach(short = dat_all$Habitat, long = hbt) %do% {
  dat_all$Habitat[dat_all$Habitat == short] <- long
}

dat_all <- dat_all[, c(1, 2, 6, 7, 3:5)]


### create table -----

stargazer(dat_all, summary = FALSE, rownames = FALSE, digits = 2)
