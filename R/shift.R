#' shifting operator
#'
#' @param im cimg object
#' @param s numeric \code{p} by \code{2} matrix containing the registration
#'   sameters
#'
#' @return shifted cimg object
#'
#' @export
#' @examples
#' shift(cameraman, c(1,1))
#' shift(cameraman, cbind(c(1,1),c(-0.5,0.5)))
shift <- function(im, s){
  if(is.null(dim(s)))
    s <- rbind(s)
  if(!is.imlist(im))
    im <- imsplit(im, 'z')
  if(length(im) == 1 & nrow(s) > 1)
    im <- as.imlist(replicate(nrow(s), im))
  if(length(im) > 1 & nrow(s) == 1)
    s <- t(replicate(length(im), s, simplify = TRUE))

  im <- lapply(1:nrow(s), function(i){
    if(any(s[i, ] %% 1 != 0))
      return(shift.subpx(im[[i]], s[i,1], s[i,2]))
    else
      return(imager::imshift(im[[i]], s[i,1], s[i,2], boundary_conditions = 1))
    })

  im <- imappend(im, 'z')
  return(im)
}

shift.subpx <- function(im, u, v){
  map <- function(x, y){
    return(list(x = u, y = v))
  }
  im <- imager::imwarp(im, map = map,
                       coordinates = "relative",
                       direction="backward",
                       interpolation = "cubic",
                       boundary = 'neumann')
  return(im)
}
