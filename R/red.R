#' REgularization by Denoising
#'
#' @param y cimg object with the low resolution frame(s)
#' @param mu numeric indicating the step size
#' @param lambda,sigma numeric indicating the regularization parameters
#' @param functional character with the optimization task or function with the functional to be used
#' @param engine character indicating the denoised engine or function with the denoiser engine to be used
#' @param niter numeric indicating the maximum number of iterations
#' @param args arguments to be passed implicity to H HT and f
#'
#' @export
#' @examples
#'
#' im <- lenna
#' y <- degrade(im, noise = 0.05)
#' x <- RED(y, sigma = 1, lambda = 5, mu = 0.1, 'DN', niter = 10)
#' par(mfrow = c(1,2), mar = c(0,0,2,0)+0.1)
#' plot(y, interp = FALSE, axes = FALSE, main = 'Degraded im')
#' mtext(paste(round(PSNR(im, y),2), 'dB'), side = 1, line = -2)
#' plot(x, interp = FALSE, axes = FALSE, main = 'Restored im')
#' mtext(paste(round(PSNR(im, x),2), 'dB'), side = 1, line = -2)
#'
#' im <- cameraman
#' y <- degrade(im, blur = 9)
#' x <- RED(y, sigma = 0.5, lambda = 4, mu = 0.2, 'DB', niter = 30)
#' par(mfrow = c(1,2), mar = c(0,0,2,0)+0.1)
#' plot(y, interp = FALSE, axes = FALSE, main = 'Degraded image')
#' mtext(paste(round(PSNR(im, y),2), 'dB'), side = 1, line = -2)
#' plot(x, interp = FALSE, axes = FALSE, main = 'Restored image')
#' mtext(paste(round(PSNR(im, x),2), 'dB'), side = 1, line = -2)
#'
#' im <- cameraman
#' L = 2
#' s <- cbind(c(0,1,2,-2,1,3,-1,-3,-1), c(0,-1,2,1,-2,-3,3,-2,-3))
#' y <- as.cimg(array(apply(s, 1, function(s) degrade(im, L = L, s = s, noise = 0.02)), c(dim(im)*c(1/L,1/L,nrow(s),1))))
#' y1 <- resize(imsplit(y,'z')[[1]], -100*L, -100*L, interpolation_type = 5)
#' x <- RED(y, sigma = 0.3, lambda = 2, mu = 0.1, functional = 'SR', niter = 30, args = list(scale = L, s=s))
#' par(mfrow = c(1,2), mar = c(0,0,2,0)+0.1)
#' plot(y1, interp = FALSE, axes = FALSE, main = 'Bicubic Interpolation')
#' mtext(paste(round(PSNR(im, y1),2), 'dB'), side = 1, line = -2)
#' plot(x, interp = FALSE, axes = FALSE, main = 'Super Resolved')
#' mtext(paste(round(PSNR(im, x),2), 'dB'), side = 1, line = -2)


RED <- function(y, lambda, sigma, mu = NULL, functional = 'SR', engine = 'MF', niter = 50, args = NULL){

  f <- NULL

  if (engine == 'MF')
    f$dn <- function(x) medianblur(x, n = 3, threshold = 0)
  else if(is.function(engine))
    f$dn <- engine
  else
    stop('Unsupported denoise engine')

  if (functional == 'SR'){
    f$H <- function(im){ #transform HR_FRAME to LR_FRAME
      im <- lapply(1:nrow(args$s), function(i){
        res <- im
        #res <- imager::imshift(res, par[i,1], par[i,2])
        res <- shift(res, args$s[i,])
        #res <- imager::isoblur(res, args$sigma_blur)
        #res <- imager::resize(res, ncol(im)/args$scale, nrow(im)/args$scale, interpolation_type = 2, boundary_conditions = 1)
        res <- imager::resize(res, ncol(im)/args$scale, nrow(im)/args$scale, interpolation_type = 5)
        return(res)
      })
      im <- imappend(as.imlist(im), 'z')
      return(im)
    }
    f$HT <- function(im){ #transform LR_FRAME to HR_FRAME
      im <- imsplit(im, 'z')
      im <- lapply(1:length(im), function(i){
        res <- im[[i]]
        #res <- imager::resize(res, ncol(im[[i]])*args$scale, nrow(im[[i]])*args$scale, interpolation_type = 2, boundary_conditions = 1)
        res <- imager::resize(res, ncol(res)*args$scale, nrow(res)*args$scale, interpolation_type = 1)
        #res <- imager::imsharpen(res, args$amplitude, type = 'shock', alpha = 5, sigma = args$sigma_blur)
        #res <- imager::imshift(res, -par[i,1], -par[i,2])
        res <- shift(res, -args$s[i,])
        return(res)
      })
      im <- parmed(as.imlist(im))
      return(im)
    }
    x <- f$HT(imsplit(y, 'z')[[1]])
  }

  if (functional == 'DN'){
    f$H <- function(im){
      return(im)
    }
    f$HT <- function(im){
      return(im)
    }
    x <- y
  }

  if (functional == 'DB'){
    if (is.null(args$filter)){
      args$filter <- imfill(9, 9, val = 1/9^2)
      #grid <- seq(-5,5,1)
      #h <- pnorm(grid + 0.5) - pnorm(grid - 0.5)
      #h <- expand.grid(h, h)
      #h <- apply(h, 1, prod)
      #args$filter <- cimg(array(h, c(11,11,1,1)))
      }
    #f$H <- function(im) return(fft_convolve(im, args$filter))
    f$H <- function(im) return(isoblur(im, sigma = 3))
    #f$HT <- function(im) return(fft_convolve(im, args$filter))
    f$HT <- function(im) return(isoblur(im, sigma = 3))
    x <- y
  }

  N <- prod(dim(x))

  if(is.null(mu))
    mu <- 2/(1/(sigma^2) + lambda)
  p <- dim(y)[3]

  ## stepest descent
  for(i in 1:niter){
    dif <- f$H(x) - y
    difn <- x - f$dn(x)

    grad <- (1/(sigma^2))*f$HT(dif) + lambda*difn
    loss <- (1/(2*(sigma^2)))*sum(dif^2)/p + (lambda/2)*sum(x*difn)
    #x <- x - mu*zero.outliers(grad)
    x <- x - mu*grad
    cat(loss/N,'\n')
#    plot(x, interp = F)

  }

  return(x)
}


### checking H and HT stability
if(F){

  ###
  args <- list(sigma_blur = 1.6, amplitude = 1.6)
  f <- NULL
  f$H <- function(im){ #transform HR_FRAME to LR_FRAME
    im <- isoblur(im, args$sigma_blur, gaussian = TRUE)
    return(im)
  }
  f$HT <- function(im){ #transform LR_FRAME to HR_FRAME
    im <- imsharpen(im, args$amplitude, type = 'shock', alpha = 5, sigma = args$sigma_blur)
    return(im)
  }

  im <- im0 <- lenna
  loss <- NULL

  for (i in 1:20){
    im <- f$HT(f$H(im))
    loss <- c(loss, MSE(im, im0))
  }
  par(mfrow = c(1,2))
  plot(loss)
  plot(im, axes = F)

  ####
  args <- list(scale = 4)
  f <- NULL
  f$H <- function(im){ #transform HR_FRAME to LR_FRAME
    im <- resize(im, ncol(im)/args$scale, nrow(im)/args$scale, interpolation_type = 5)
    return(im)
  }
  f$HT <- function(im){ #transform LR_FRAME to HR_FRAME
    im <- resize(im, ncol(im)*args$scale, nrow(im)*args$scale, interpolation_type = 1)
    return(im)
  }

  im <- im0 <- lenna
  loss <- NULL

  for (i in 1:20){
    im <- f$HT(f$H(im))
    loss <- c(loss, MSE(im, im0))
  }
  par(mfrow = c(1,2))
  plot(loss)
  plot(im, axes = F)

  ####
  args <- list(par = c(1,1))
  f <- NULL
  f$H <- function(im){ #transform HR_FRAME to LR_FRAME
    im <- shift(im, par = args$par)
    return(im)
  }
  f$HT <- function(im){ #transform LR_FRAME to HR_FRAME
    im <- shift(im, par = - args$par)
    return(im)
  }

  im <- im0 <- lenna
  loss <- NULL

  for (i in 1:20){
    im <- f$HT(f$H(im))
    loss <- c(loss, MSE(im, im0))
  }
  par(mfrow = c(1,2))
  plot(loss)
  plot(im, axes = F)

  ####
  args <- list(par = c(0.5,0.5))
  f <- NULL
  f$H <- function(im){ #transform HR_FRAME to LR_FRAME
    im <- shift(im, par = args$par)
    return(im)
  }
  f$HT <- function(im){ #transform LR_FRAME to HR_FRAME
    im <- shift(im, par = - args$par)
    return(im)
  }

  im <- im0 <- lenna
  loss <- NULL

  for (i in 1:20){
    im <- f$HT(f$H(im))
    loss <- c(loss, MSE(im, im0))
  }
  par(mfrow = c(1,2))
  plot(loss)
  plot(im, axes = F)

}

if(F){
  Lambda <- 10^(seq(-2,1,0.25))
  y <- degrade(lenna, sd_noise = 0.05, sd_blur = 0)
  x <- NULL
  mse <- MSE(y, lenna)
  for(lambda in Lambda){
    im <- RED(y, sigma = 1, lambda = lambda, mu = 0.1, 'DN', niter = 10)
    mse <- c(mse, MSE(im, lenna))
    x <- c(x, imsplit(im, 'z'))
  }
  plot(c(0,Lambda), mse)
  display(as.imlist(x))
}


