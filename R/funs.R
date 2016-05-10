#' Compute Root Mean Square error
#' 
rmse <- function(x, y) {
  sqrt(mean((y - x)^2, na.rm = TRUE))
}

#' Setup Moving Window
#' 
movingWindow <- function(width = 7L, height = width) {
  mat <- matrix(1, nrow = height, ncol = width)
  mat[ceiling(height/2), ceiling(width/2)] <- 0
  mat
}