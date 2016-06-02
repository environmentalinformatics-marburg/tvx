lower <- rev(sortElevation(FALSE))[1:5]
middle <- rev(sortElevation(FALSE))[6:10]
upper <- sortElevation(FALSE)[1:3]

dat_training[dat_training$habitat %in% lower, ]
dat_training[dat_training$habitat %in% middle, ]
dat_training[dat_training$habitat %in% upper, ]

### -----

prediction_stats[prediction_stats$habitat %in% lower, ]
prediction_stats[prediction_stats$habitat %in% middle, ]
prediction_stats[prediction_stats$habitat %in% upper, ]


### -----

rst_r <- raster("data/results/ndvi_lst_r.tif")
rst_p <- raster("data/results/ndvi_lst_p.tif")

mean(rst_r[rst_p[] < 0.001])
range(rst_r[rst_p[] < 0.001])
