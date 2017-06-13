#' Convolution of two images via FFT
#'
#' @param im,filter \code{cimg} objects
#' @param deconvolution logical indicating if the deconvolution should be performed
#'
#' @export
#' @examples
#' im <- lenna
#' filter <- imfill(9,9,val = 1)
#' blurred.im <- fft_convolve(im, filter)
#' deblurred.im <- fft_convolve(blurred.im, filter, deconvolution = TRUE)
#' par(mfrow = c(1,3), mar = c(0,0,1,0)+0.1)
#' plot(im, axes = FALSE, interp = FALSE, main = 'Original Lenna')
#' plot(blurred.im, axes = FALSE, interp = FALSE, main = 'Blurred Lenna')
#' plot(deblurred.im, axes = FALSE, interp = FALSE, main = 'deBlurred Lenna')
#' PSNR(im, blurred.im)
#' PSNR(im, deblurred.im)
fft_convolve <- function(im, filter, deconvolution = FALSE){

  dim <- dim(im)
  if(dim[1] != dim[2])
    stop("'im' must be square")
  dim <- dim[1]
  dfi <- dim(filter)
  if(dfi[1] != dfi[2])
    stop("'filter' must be square")
  dfi <- dfi[1]

  filter <- filter/sum(filter)
  if(dfi %% 2 == 1)
    filter <- imager::pad(filter, 1, 'xy', 1)
  filter <- imager::pad(filter, dim - dfi - 1, 'xy')

#  filter <- imager::pad(filter, 2*dfi, 'xy')
#  im <- imager::pad(im, 2*dfi, 'xy')

  im <- stats ::fft(im)
  filter <- stats ::fft(filter)

  if(deconvolution){
    res <- im / filter #https://en.wikipedia.org/wiki/Wiener_deconvolution ??? flip kernel???
  }
  else
    res <- im * filter

  res <- stats ::fft(res, inverse = TRUE) / length(res)
  res <- Re(res)
  res <- imappend(imsplit(res, 'x', 2)[2:1], 'x')
  res <- imappend(imsplit(res, 'y', 2)[2:1], 'y')

#  res[px.borders(res, dfi)] <- 0 # employ the option for dirishlet borders
#  res <- autocrop(res)

  return(res)
}
