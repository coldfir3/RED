#' Registration parameter estimation
#'
#' @param src,tar cimg objects
#' @param method character indicating the method to be used
#' @param par0 numeric vector for the initial guess for the registration parameters
#' @param verbosity Numeric indicating the level of verbosity is displayed
#' @param ... parameters to be passed to the optimization algorithm
#'
#' @return the registration parameters, usualy a 2d vector.
#'
#' @export
#' @examples
#'
#' src <- cameraman
#'
#' tar <- shift(cameraman, c(5,-15))
#' round(s <- register(src, tar, method = 'coarse', steps = 4), 4)
#'
#' tar <- shift(cameraman, c(-1.155, 3.231))
#' round(s <- register(src, tar, method = 'taylor', tol = 1e-4), 4)
#'
#' tar <- transform(cameraman, c(c(-1.155, 1.231, 0.121)))
#' round(s <- register(src, tar, method = 'taylor3', tol = 1e-4, maxiter = 100), 4)
register <- function(src, tar, method = 'taylor', par0 = c(0,0,0), verbosity = 2, ...){

  if (method == 'taylor')
    par <- reg.taylor(src, tar, par0 = par0[1:2], verbosity = verbosity, ...)
  else if (method == 'taylor3')
    par <- reg.taylor3(src, tar, par0 = par0[1:3], verbosity = verbosity, ...)
  else if (method == 'coarse')
    par <- reg.NSS(src, tar, par0 = par0[1:2], verbosity = verbosity, ...)
  else
    stop("method not available.")

  return(par)
}

reg.taylor <- function(src, tar, par0 = c(0,0), verbosity = 2, ...){

  args <- list(...)
  if(is.null(args$maxiter))
    maxiter = 50
  else
    maxiter = args$maxiter
  if(is.null(args$tol))
    tol = 1e-3
  else
    tol = args$tol

  par.acc <- par <- par0
  src.ori <- src
  tar.ori <- tar

  if(verbosity > 1)
    cat("\nDetermining registration parameters by Taylor series:\n")
  if(verbosity > 1)
    cat("iteration number:", 0,"\t| Par:", round(par.acc,3), '\tMSE:', NA, '\n')
  for (i in 1:maxiter){

    src <- shift(src.ori, par.acc)

    dif <- tar - src
    src_g <- imager::imgradient(src, "xy")

    A11 <- sum(src_g$x^2)
    A12 <- A21 <- sum(src_g$x * src_g$y)
    A22 <- sum(src_g$y^2)

    A <- matrix(c(A11, A21, A12, A22), nrow = 2)
    b <- c(sum(src_g$x * dif), sum(src_g$y * dif))

    par <- - solve(A, b)
    par.acc <- par.acc + par

    err <- MSE(dif)

    if (i %in% round(2^seq(0, log2(maxiter), length.out = 6)) & verbosity > 1)
      cat("iteration number:", i,"\t| Par:", round(par.acc,3), '\tMSE:', round(err, 3), '\n')

    if (max(abs(par)) < tol)
      break
  }
  if(verbosity > 0)
    cat("iteration number:", i,"\t| Par:", round(par.acc,3), '\tMSE:', round(err, 3), '(FINAL)\n')
  return(par.acc)
}

reg.taylor3 <- function(src, tar, par0 = c(0,0,0), verbosity = 2, ...){

  args <- list(...)
  if(is.null(args$maxiter))
    maxiter = 50
  else
    maxiter = args$maxiter
  if(is.null(args$tol))
    tol = 1e-3
  else
    tol = args$tol

  par.acc <- par <- par0
  src.ori <- src
  tar.ori <- tar

  if(verbosity > 1)
    cat("\nDetermining registration parameters by Taylor series:\n")
  if(verbosity > 1)
    cat("iteration number:", 0,"\t| Par:", round(par.acc,3), '\tMSE:', NA, '\n')
  for (i in 1:maxiter){

    src <- transform(src.ori, par.acc)

    dif <- tar - src
    g <- imager::imgradient(src, "xy")

    xc <- (ncol(tar) + 1)/2
    yc <- (nrow(tar) + 1)/2

    xy <- as.data.frame(tar)
    x <- xy$x - xc
    y <- xy$y - yc
    a <- x*g$y - y*g$x

    A11 <- sum(g$x^2)
    A12 <- sum(g$x * g$y)
    #    A13 <- sum(a*g$x)
    #    A23 <- sum(a*g$y)
    A22 <- sum(g$y^2)
    #    A33 <- sum(a^2)

    #    A <- matrix(c(A11, A12, A13, A12, A22, A23, A13, A23, A33), nrow = 3)
    #    b <- c(sum(g$x*dif), sum(g$y*dif), sum(a*dif))
    A <- matrix(c(A11, A12, A12, A22), nrow = 2)
    b <- c(sum(g$x*dif), sum(g$y*dif))

    par <- - solve(A, b)
    par <- c(par, -sum(a*dif)/sum(a^2) * 180/pi)
    par.acc <- par.acc + par

    err <- MSE(dif)

    if (i %in% round(2^seq(0, log2(maxiter), length.out = 6)) & verbosity > 1)
      cat("iteration number:", i,"\t| Par:", round(par.acc,3), '\tMSE:', round(err, 3), '\n')

    if (max(abs(par)) < tol)
      break
  }
  if(verbosity > 0)
    cat("iteration number:", i,"\t| Par:", round(par.acc,3), '\tMSE:', round(err, 3), '(FINAL)\n')
  return(par.acc)
}

reg.NSS <- function(src, tar, par0 = c(0,0), ...){ ## melhorar o algoritmo de otimização

  args <- list(...)
  if(is.null(args$steps))
    steps = 3
  else
    steps = args$steps

  xc <- par0[1]
  yc <- par0[2]
  S <- 2^steps
  hist <- NULL
  while(S >= 0.5){
    s <-  unique(data.frame(
      dx=c(seq(-round(S+1e-3), round(S+1e-3), 1) + xc, rep(xc, 2*round(S+1e-3) + 1)),
      dy=c(rep(yc, 2*round(S+1e-3) + 1), seq(-round(S+1e-3), round(S+1e-3), 1) + yc))
    )
    if(!is.null(hist))
      s <- unique(rbind(hist[,1:2],s))[-(1:nrow(hist)),]
    hist <- rbind(hist, cbind(s, mse = apply(s, 1, function(s) MSE(shift(src, s), tar))))
    xc <- hist[which.min(hist$mse),1]
    yc <- hist[which.min(hist$mse),2]
    S <- S/2
    #plot(hist[,1:2])
    #points(s, col = 'green')
    #points(xc,yc, col = 'red')
  }
  par <- c(xc, yc)
  return(par)
}
