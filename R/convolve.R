### deconvolução!!!! PORRAAA FILHA DA PUTAAA FOIII!!!!!
if(F){
  filter <- imfill(9,9,val = 1)
  filter <- filter/sum(filter)
  im <- im0 <- lenna

  ffilter <- pad(filter, 1, 'xy', 1)
  ffilter <-pad(ffilter, nrow(im)-10, 'xy', 0)
  fft.filter <- FFT(ffilter)
  fft.filter <- fft.filter$real + fft.filter$imag*1i
  fft.im <- FFT(im)
  fft.im <- fft.im$real + fft.im$imag*1i
  fft.conv <- fft.filter * fft.im
  conv <- FFT(Re(fft.conv), Im(fft.conv), inverse = TRUE)$real
  conv <- imappend(imsplit(imappend(imsplit(conv, 'x', 2)[2:1], 'x'), 'y', 2)[2:1], 'y')
  plot(conv, interp = FALSE)
  plot(convolve(im, filter), interp = FALSE)

  im <- conv
  fft.im <- FFT(im)
  fft.im <- fft.im$real + fft.im$imag*1i
  fft.deconv <- fft.im / fft.filter
  deconv <- FFT(Re(fft.deconv), Im(fft.deconv), inverse = TRUE)$real
  deconv <- imappend(imsplit(imappend(imsplit(deconv, 'x', 2)[2:1], 'x'), 'y', 2)[2:1], 'y')
  plot(deconv, interp = FALSE)
  plot(im0, interp = FALSE)
}


#' Convolution of two images via FFT
#'
#' @param im,filter cimg objects
#'
#' @export
#' @examples
#' im <- lenna
#' filter <- imfill(9,9,val = 1)
#' blurred.im <- fft_conv(im, filter)
#' deblurred.im <- fft_conv(blurred.im, filter, deconvolution = TRUE)
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

  im <- fft(im)
  filter <- fft(filter)

  if(deconvolution){
    res <- im / filter #https://en.wikipedia.org/wiki/Wiener_deconvolution ??? flip kernel???
  }
  else
    res <- im * filter
  res <- fft(res, inverse = TRUE) / length(res)
  res <- Re(res)
  res <- imappend(imsplit(res, 'x', 2)[2:1], 'x')
  res <- imappend(imsplit(res, 'y', 2)[2:1], 'y')

#  res[px.borders(res, dfi)] <- 0 # employ the option for dirishlet borders
#  res <- autocrop(res)

  return(res)
}
