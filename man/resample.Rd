% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resample.R
\name{resample}
\alias{resample}
\title{Resampling of an image}
\usage{
resample(im, L = 1, L1 = L, L2 = L)
}
\arguments{
\item{im}{\code{cimg} object}

\item{L}{numeric indicating the overall scale change. This parameter will be
override by L1 or L2}

\item{L1, L2}{numeric indicating the directional scale change}
}
\value{
A resampled \code{cimg} object
}
\description{
Resampling of an image
}
\examples{
im <- lenna
par(mfrow = c(1,2), mar = rep(0,4)+0.1)
plot(im, axes = FALSE, interp = FALSE)
plot(resample(im, 1/4), axes = FALSE, interp = FALSE)
}
