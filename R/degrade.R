#' Degradation of an image
#'
#' This function degrades a high resolution image into a low resolution
#' image. The degradation folows this: bla bla bla
#'
#' @param z a cimg object containing the high resolution image
#' @inheritParams resample
#' @inheritParams shift
#' @param noise numeric indicating the standard deviation of the noise or an cimg object that will be added to the resampled z
#' @param blur numeric indicating the blur range (for uniform blur) or an cimg object with the blur kernel to be convolved with z
#'  if nothing is provided an default kernel will be used.
#'
#' @return A degraded cimg object
#'
#' @export
#' @examples
#' degraded.lenna <- degrade(lenna, L = 4, noise = 0.05, blur = 3)
#' par(mfrow = c(1,2), mar = c(0,0,1,0)+0.1)
#' plot(lenna, axes = FALSE, interp = FALSE, main = 'Original Lenna')
#' plot(degraded.lenna, axes = FALSE, interp = FALSE, main = 'Degraded Lenna')
degrade <- function(z, L = 1, s = cbind(0,0), noise = 0, blur = 1, L1 = L, L2 = L){

  p <- nrow(s)
  N <- dim(z)
  M <- round(c(N[1]/L1, N[2]/L2, p, 1))

  if (is.cimg(noise))
    noise <- noise
  else if (noise == 0)
    noise <- noise
  else if(noise > 0)
    noise <- imager::as.cimg(array(stats::rnorm(prod(M), mean = 0, sd = noise), M))
  else
    stop("noise must be non-negative scalar or cimg object")

  if (is.cimg(blur))
    blur <- blur
  else if (length(blur) == 1)
    blur <- imager::imfill(blur, blur, 1, 1/blur^2)
  else if (is.null(blur)) # change condition to a name
    blur <- imager::imfill(5, 5, 1, 1/5^2) #change to PSF #incoherent psf finite detector

  y <- z
  y <- shift(y, s)
  y <- convolve(y, blur, dirichlet = TRUE)
  y <- resample(y, L1 = 1/L1, L2 = 1/L2)
  y <- y + noise

  return(y)
}
