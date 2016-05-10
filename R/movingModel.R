movingModel <- function(lst, ndvi, directions = movingWindow(), 
                        required = 2/3, rasterize = TRUE, ...) {
  
  ## sample size
  n <- prod(dim(directions))

  ## convert layers to matrices
  mat_x <- t(as.matrix(ndvi))
  mat_y <- t(as.matrix(lst))
  
  ## loop over cells
  metrics <- lapply(1:ncell(lst), function(i) {
    
    # identify adjacent cells
    id <- adjacent(lst, i, directions, include = TRUE, pairs = FALSE, ...)
    
    # stop if number of adjacent cells is below threshold; this typically
    # affects edge pixels only
    if (length(id) < required * n) {
      return(data.frame(cell = i, slope = NA, intercept = NA, r = NA))
    }
    
    # retrieve cell values
    val_x <- mat_x[id] 
    val_y <- mat_y[id]
    
    dat <- data.frame(x = val_x, y = val_y)
    
    # proceed if number of missing values does not exceed threshold
    if (sum(complete.cases(dat)) >= required * n) {
      dat <- dat[complete.cases(dat), ]
      
      # calculate slope, intercept and correlation coefficient
      r <- cor(dat$x, dat$y)
      mod <- lm(y ~ x, data = dat)
      
      data.frame(cell = i, slope = coef(mod)[2], intercept = coef(mod)[1], r = r)
      
    } else {
      data.frame(cell = i, slope = NA, intercept = NA, r = NA)
    }
  })
  
  metrics <- do.call("rbind", metrics)
  
  ## rasterize data (optional)
  if (rasterize) {
    metrics <- lapply(c("slope", "intercept", "r"), function(i) {
      setValues(lst, metrics[, i])
    })
    
    names(metrics) <- c("slope", "intercept", "r")
  }
  
  return(metrics)
}
