#' REgularization by Denoising
#'
#' @param y cimg object with the observed frame(s)
#' @param x0 initial guess for the output image, if NULL an educated guess will
#'   be used. If a custom functional is provided this cant be NULL
#' @param step numeric indicating the step size (if NULL an optimal step size
#'   will be used)
#' @param lambda,sigma numeric indicating the regularization parameters
#' @param tol numeric indicating the stoping criteria. The algorithm will stop
#'   when step < tol. Default = 0.001
#' @param functional character with the optimization task or function with the
#'   functional to be used
#' @param engine character indicating the denoised engine or function with the
#'   denoiser engine to be used
#' @param niter numeric indicating the maxistepm number of iterations
#' @param args arguments to be passed implicity to H HT and f
#'
#' @export
#' @examples
#'
#' im <- lenna
#' y <- degrade(im, noise = 0.05)
#' x <- RED(y, sigma = 1, lambda = 5, functional = 'DN', niter = 50)
#' par(mfrow = c(1,2), mar = c(0,0,2,0)+0.1)
#' plot(y, interp = FALSE, axes = FALSE, main = 'Degraded im')
#' mtext(paste(round(PSNR(im, y),2), 'dB'), side = 1, line = -2)
#' plot(x, interp = FALSE, axes = FALSE, main = 'Restored im')
#' mtext(paste(round(PSNR(im, x),2), 'dB'), side = 1, line = -2)
#'
#' \dontrun{
#' im <- cameraman
#' y <- degrade(im, blur = 5)
#' y<- isoblur(im, 3, gaussian = TRUE)
#' x <- RED(y, sigma = 1, lambda = 4, functional = 'DB', niter = 1500)
#' par(mfrow = c(1,2), mar = c(0,0,2,0)+0.1)
#' plot(y, interp = FALSE, axes = FALSE, main = 'Degraded image')
#' mtext(paste(round(PSNR(im, y),2), 'dB'), side = 1, line = -2)
#' plot(x, interp = FALSE, axes = FALSE, main = 'Restored image')
#' mtext(paste(round(PSNR(im, x),2), 'dB'), side = 1, line = -2)
#'
#' im <- cameraman
#' L = 2
#' s <- cbind(c(0,1,2,-2,1,3,-1,-3,-1), c(0,-1,2,1,-2,-3,3,-2,-3))
#' y <- degrade(im, L = L, s = s, noise = 0.05)
#' xref <- resize(imsplit(y,'z')[[1]], -100*L, -100*L, interpolation_type = 5)
#' x <- RED(y, sigma = 1, lambda = 5, functional = 'SR', niter = 50, args = list(scale = L, s=s))
#' par(mfrow = c(1,2), mar = c(0,0,2,0)+0.1)
#' plot(xref, interp = FALSE, axes = FALSE, main = 'Bicubic Interpolation')
#' mtext(paste(round(PSNR(im, xref),2), 'dB'), side = 1, line = -2)
#' plot(x, interp = FALSE, axes = FALSE, main = 'Super Resolved')
#' mtext(paste(round(PSNR(im, x),2), 'dB'), side = 1, line = -2)
#'
#' im0 <- 0.2*pad(cameraman, 256, 'xy')
#' im1 <- lenna
#' im2 <- im1 - im0
#' y1 <- degrade(im1, noise = 0.05)
#' y2 <- degrade(im2, noise = 0.05)
#' y0 <- y1 - y2
#' x0 <- RED(y0, sigma = 1, lambda = 50, functional = 'DN', niter = 100)
#'
#' par(mfrow = c(1,2), mar = c(0,0,2,0)+0.1)
#' plot(y0, interp = FALSE, axes = FALSE, main = 'naive')
#' mtext(paste(round(PSNR(im0, y0),2), 'dB'), side = 1, line = -2)
#' plot(x0, interp = FALSE, axes = FALSE, main = 'proposed')
#' mtext(paste(round(PSNR(im0, x0),2), 'dB'), side = 1, line = -2)
#' }
RED <- function(y, x0 = NULL, lambda = 1, sigma = 1, functional = 'SR', engine = 'MF', niter = 50, step = NULL, tol = 1e-3, args = NULL){

  f <- NULL

  if (engine == 'MF'){
    if (is.null(args$n)){
      if(functional == 'SR')
        args$n <- args$scale+1
      else
        args$n <- 3
      }
      f$dn <- function(x) medianblur(x, n = args$n, threshold = 0)
    }
  else if(is.function(engine))
    f$dn <- engine
  else
    stop('Unsupported denoise engine')

  if (functional == 'SR'){
    if(is.null(args$s))
      stop('Relocation parameters must be provided via "args = list(s = )"')
    if(is.null(args$scale))
      stop('Scale parameter must be provided via "args = list(scale = )"')
    if(nrow(args$s) != dim(y)[3])
      stop('Number of relocation parameters are not compatible with number of frames')
    f$H <- function(im){
      im <- lapply(1:nrow(args$s), function(i){
        res <- im
        res <- shift(res, args$s[i,])
        res <- resample(res, 1/args$scale)
        return(res)
      })
      im <- imappend(as.imlist(im), 'z')
      return(im)
    }
    f$HT <- function(im){
      im <- imsplit(im, 'z')
      im <- lapply(1:length(im), function(i){
        res <- im[[i]]
        res <- resample(res, args$scale)
        res <- shift(res, -args$s[i,])
        return(res)
      })
      im <- parmed(as.imlist(im))
      return(im)
    }
    if(is.null(x0))
      x <- f$HT(imsplit(y, 'z')[[1]])
    else
      x <- x0
  }

  if (functional == 'DN'){
    f$H <- function(im){
      return(im)
    }
    f$HT <- function(im){
      return(im)
    }
    if(is.null(x0))
      x <- y
    else
      x <- x0
  }

  if (functional == 'DB'){
    f$H <- function(im) return(isoblur(im, sigma = 3, gaussian = TRUE))
    f$HT <- function(im) return(isoblur(im, sigma = 3, gaussian = TRUE))
    if(is.null(x0))
      x <- y
    else
      x <- x0
  }

  if(is.null(step)){
    step <- 2/(1/(sigma^2) + lambda)
    step_iter <- TRUE
    }
  else
    step_iter <- FALSE

  dif <- f$H(x) - y
  difn <- x - f$dn(x)
  mse <- (1/(2*(sigma^2)))*sum(dif^2)
  penalty <- (lambda/2)*sum(x*difn)
  loss <- mse + penalty

  ## stepest descent
  for(n in 1:niter){

    grad <- (1/(sigma^2))*f$HT(dif) + lambda*difn
    x <- x - step*grad
    dif <- f$H(x) - y
    difn <- x - f$dn(x)
    mse <- c((1/(2*(sigma^2)))*sum(dif^2), mse)
    penalty <- c((lambda/2)*sum(x*difn), penalty)
    loss <- c(mse[1] + penalty[1], loss)

    if(loss[1] > loss[2] & step_iter){
      step <- step/2
      x <- x + step*grad
    }
    else if(step_iter)
      step <- 1.1*step
    if(step < tol)
      break

    cat("\niteration:", n, "of", niter, "(max) | step size:", round(step,3),
        "| Loss:", round(loss[1]/utils::tail(loss,1),3),
        'MSE:', round(mse[1]/loss[1], 3),
        'penalty:', round(penalty[1]/loss[1], 3))
  }
  return(x)
}

if(F){ #gaussian PSF

  #grid <- seq(-5,5,1)
  #h <- pnorm(grid + 0.5) - pnorm(grid - 0.5)
  #h <- expand.grid(h, h)
  #h <- apply(h, 1, prod)
  #args$filter <- cimg(array(h, c(11,11,1,1)))
}


