#' Transform an image
#'
#' @param im \code{cimg} object
#' @param s numeric \code{1} by \code{3} vector containing the registration
#'   parameters
#'
#' @return shifted \code{cimg} object
#'
#' @export
#' @examples
#' shift(cameraman, c(1,1))
#' shift(cameraman, cbind(c(1,1),c(-0.5,0.5)))
transform <- function(im, s){ # change s to a matrix (see shift function)

  theta <- s[3] * pi/180
  u <- s[1]
  v <- s[2]

  xC <- (ncol(im)+1)/2
  yC <- (nrow(im)+1)/2
  sn <- sin(theta)
  cs <- cos(theta)
  map <- function(x, y){
    xl <- x - xC
    yl <- y - xC
    x <-  cs * xl + sn * yl - u + xC
    y <- -sn * xl + cs * yl - v + xC
    return(list(x = x, y = y))
  }
  im <- imager::imwarp(im, map = map, direction="backward", interpolation = "cubic")

  return(im)
}
