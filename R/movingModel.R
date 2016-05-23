movingModel <- function(lst, ndvi, dem, directions = movingWindow(), 
                        required = 2/3, rasterize = TRUE, ...) {
  
  ## sample size
  n <- prod(dim(directions))
  
  ## convert lst, ndvi layers to matrices
  mat_x <- t(as.matrix(ndvi))
  mat_y <- t(as.matrix(lst))
  
  ### without dem -----
  if (missing(dem)) {
    
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
    
    
  ### including dem -----    
  } else {
    
    ## convert dem to matrix
    mat_z <- t(as.matrix(dem))
    
    ## loop over cells
    metrics <- lapply(1:ncell(lst), function(i) {
      
      # identify adjacent cells
      id <- adjacent(lst, i, directions, include = TRUE, pairs = FALSE, ...)
      
      # stop if number of adjacent cells is below threshold; this typically
      # affects edge pixels only
      if (length(id) < required * n) {
        return(data.frame(cell = i, slope1 = NA, slope2 = NA, intercept = NA, 
                          r1 = NA, r2 = NA))
      }
      
      # retrieve cell values
      val_x <- mat_x[id] 
      val_y <- mat_y[id]
      val_z <- mat_z[id]
      
      dat <- data.frame(x = val_x, y = val_y, z = val_z)
      
      # proceed if number of missing values does not exceed threshold
      if (sum(complete.cases(dat)) >= required * n) {
        dat <- dat[complete.cases(dat), ]
        
        # calculate slope, intercept and correlation coefficient
        r <- cor(dat)
        mod <- lm(y ~ x + z, data = dat)
        
        data.frame(cell = i, slope1 = coef(mod)[2], slope2 = coef(mod)[3], 
                   intercept = coef(mod)[1], r1 = r[1, 2], r2 = r[2, 3])
        
      } else {
        data.frame(cell = i, slope1 = NA, slope2 = NA, intercept = NA, 
                   r1 = NA, r2 = NA)
      }
    })
    
    metrics <- do.call("rbind", metrics)
    
    ## rasterize data (optional)
    if (rasterize) {
      metrics <- lapply(c("slope1", "slope2", "intercept", "r1", "r2"), 
                        function(i) {
                          setValues(lst, metrics[, i])
                        })
      
      names(metrics) <- c("slope1", "slope2", "intercept", "r1", "r2")
    }
  }
  
  return(metrics)
}
