### coffee -----

fls_cof3 <- list.files("data/station_data", pattern = "^plot_cof3.*.csv", 
                       full.names = TRUE)

lst_cof3 <- lapply(fls_cof3, function(i) {
  read.csv(i)[, c("datetime", "Ta_200")]
})

dat_cof3 <- Reduce(function(...) merge(..., by = 1, all = TRUE), lst_cof3)

dat_cof3_mrg <- data.frame("plotID" = "cof3", 
                           "datetime" = dat_cof3$datetime, 
                           "Ta_200" = rowMeans(dat_cof3[, 2:4], na.rm = TRUE))

write.csv(dat_cof3_mrg, "data/station_data/plots/cof3.csv", 
          row.names = FALSE, quote = FALSE)

### grassland -----

fls_gra1 <- "data/station_data/plot_gra1_51021020181.csv"
dat_gra1 <- read.csv(fls_gra1)

write.csv(data.frame(plotID = "gra1", dat_gra1[, c("datetime", "Ta_200")]), 
          "data/station_data/plots/gra1.csv", row.names = FALSE, quote = FALSE)