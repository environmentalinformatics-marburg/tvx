#' Setup Moving Window
#' 
movingWindow <- function(width = 7L, height = width) {
  mat <- matrix(1, nrow = height, ncol = width)
  mat[ceiling(height/2), ceiling(width/2)] <- 0
  mat
}
