crop_borders <- function(im, N, center){
  x <- y <- NULL #this line is only here to prevent a "NOTE" when builing the package
  center <- round(center)
  if(center[1] - N/2 < 0){
    lx <- 0
    hx <- N
  }
  else if(center[1] + N/2 - 1 > ncol(im)){
    lx <- ncol(im) - N
    hx <- ncol(im) - 1
  }
  else{
    lx <- center[1] - N/2
    hx <- center[1] + N/2 - 1
  }
  if(center[2] - N/2 < 0){
    ly <- 0
    hy <- N
  }
  else if(center[2] + N/2 - 1 > nrow(im)){
    ly <- nrow(im) - N
    hy <- nrow(im) - 1
  }
  else{
    ly <- center[2] - N/2
    hy <- center[2] + N/2 - 1
  }

  im <- imager::imsub(im, x %inr% c(lx, hx), y %inr% c(ly, hy))
  return(im)
}
