#' Resampling of an image
#'
#' @param im \code{cimg} object
#' @param L numeric indicating the overall scale change. This param will be
#'   override by L1 or L2
#' @param L1,L2 numeric indicating the directional scale change
#'
#' @return A resampled \code{cimg} object
#'
#' @export
#' @examples
#' im <- lenna
#' par(mfrow = c(1,2), mar = rep(0,4)+0.1)
#' plot(im, axes = FALSE, interp = FALSE)
#' plot(resample(im, 1/4), axes = FALSE, interp = FALSE)
resample <- function(im, L = 1, L1 = L, L2 = L){
  N <- dim(im)
  M <- round(c(N[1]*L1, N[2]*L2, N[3], N[4]))
  if(L1 > 1 & L2 > 1)
    im <- imager::resize(im, M[1], M[2], interpolation_type = 1)
  else if(L1 < 1 & L2 < 1)
    im <- imager::resize(im, M[1], M[2], interpolation_type = 5)
  else if(L1 == 1 & L2 == 1)
    im <- im
  else if(L1 > 1)
    im <- imager::resize(im, M[1], N[2], interpolation_type = 1)
  else if(L1 < 1)
    im <- imager::resize(im, M[1], N[2], interpolation_type = 5)
  else if(L2 > 1)
    im <- imager::resize(im, N[1], M[2], interpolation_type = 1)
  else if(L2 < 1)
    im <- imager::resize(im, N[1], M[2], interpolation_type = 5)

  return(im)
}
