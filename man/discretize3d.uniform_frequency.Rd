% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discretization.plugin.R
\name{discretize3d.uniform_frequency}
\alias{discretize3d.uniform_frequency}
\title{Discretize three sets of continuous data into 3-dimensional bins by uniform frequency}
\usage{
discretize3d.uniform_frequency(x, y, z)
}
\arguments{
\item{x}{a numeric vector of the first random variable.}

\item{y}{a numeric vector of the second random variable.}

\item{z}{a numeric vector of the third random variable.}
}
\value{
discretize3d.uniform_frequency returns a 3-dimensional count table.
}
\description{
discretize3d.uniform_frequency assigns the observations of three continuous random variables to bins according to the "uniform frequency" method,
and returns a corresponding 3-dimensional count table.
}
\details{
Uniform frequency-based method ("uniform_frequency") divides the continuous data into N bins with (approximate) equal count number.
The number of bins N is initialized into a round-off value according to the square root of the data size.
}
\examples{
# three numeric vectors corresponding to three continuous random variables
x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)

# corresponding joint count table estimated by "uniform frequency" algorithm
discretize3d.uniform_frequency(x,y,z)
}
