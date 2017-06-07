#' Error measurements of images
#'
#' This function calculates error between two images
#'
#' @param x,y \code{cimg} objects
#'
#' @examples
#' degraded.lenna <- degrade(lenna, noise = 0.05)
#' MSE(lenna, degraded.lenna)
#' MAE(lenna, degraded.lenna)
#' #alternatively it can be done like:
#' MSE(lenna - degraded.lenna)
#' MAE(lenna - degraded.lenna)
#' @name error
NULL

#' @export
#' @describeIn error Mean Squared Error
MSE <- function(x, y = NULL){

  if(is.null(y))
    return(mean(x^2))
  else
    return(mean((x - y)^2))
}

#' @export
#' @describeIn error Mean Absolute Error
MAE <- function(x, y = NULL){

  if(is.null(y))
    return(mean(abs(x)))
  else
    return(mean(abs(x - y)))
}
