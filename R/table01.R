### environmental stuff -----

## clear workspace
rm(list = ls(all = TRUE))

## load packages
lib <- c("rgdal", "dplyr", "foreach", "stargazer")
Orcs::loadPkgs(lib)

## load functions
source("R/funs.R")


### ndvi_max -----


dat_all <- readRDS("data/results/stats_intra_7by7_terra.rds")[, 1:2]
dat_all$NDVImax <- round(dat_all$NDVImax, 2)
dat_all$"Training plots" <- c("cof1,3", "fed3,4", "fer0", "flm1,4", "foc1,2", 
                              "fod3,5", "fpd3,5", "fpo4,5", "gra4,5", 
                              "hel4,5", "hom4,5", "mai3,4", "sav1,3")

dat_all <- dat_all[, c(1, 3, 2)]
rownames(dat_all) <- NULL
names(dat_all)[1] <- "Habitat"


### elevations -----

dat_ele <- sortElevation()
dat_all <- merge(dat_ele, dat_all, by = "Habitat", all = TRUE, sort = FALSE)


### disturbance status -----

## disturbed habitats
dst <- c("fed", "fpd", "fod", "gra", "hom", "cof", "mai")

dat_all$Disturbance <- "No"
dat_all$Disturbance[dat_all$Habitat %in% dst] <- "Yes"
dat_all <- dat_all[, c("Habitat", "Disturbance", "Elevation", "Training plots", 
                       "NDVImax")]


### long habitat names -----

hbt <- c("Helichrysum", "Erica forest disturbed", "Erica forest", 
         "Podocarpus forest disturbed", "Podocarpus forest", 
         "Ocotea forest disturbed", "Ocotea forest", 
         "Lower montane forest", "Grassland", "Chagga homegarden", 
         "Coffee plantation", "Savanna", "Maize field")

jnk <- foreach(short = dat_all$Habitat, long = hbt) %do% {
  dat_all$Habitat[dat_all$Habitat == short] <- long
}


### create table -----

stargazer(dat_all, summary = FALSE, rownames = FALSE, digits = 2)
