#' Discretize a set of continuous data into 1-dimensional bins by "uniform width" method
#'
#' discretize1d.uniform_width assigns the observations of continuous random variables to bins according to the "uniform width" method,
#' and returns a corresponding count table.
#'
#' Uniform width-based method ("uniform_width") divides the continuous data into N bins with equal width.
#' The number of bins N is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of a random variable.
#'
#' @return discretize1d.uniform_width returns a count table.
#' @export
#'
#' @examples
#'
#' # a numeric vector corresponding to a continuous random variable
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' discretize1d.uniform_width(x)
discretize1d.uniform_width = function(x)
{
  numBins = floor( sqrt(length(x)) )
  r=range(x)

  b = seq(from=r[1], to=r[2], length.out=numBins+1 )
  X = table( cut(x, breaks=b , include.lowest=TRUE) )

  return(X)
}


#' Discretize two sets of continuous data into 2-dimensional bins by "uniform width" method
#'
#' discretize2d.uniform_width assigns the observations of two continuous random variables to bins according to the "uniform width" method,
#' and returns a corresponding 2-dimensional count table.
#'
#' Uniform width-based method ("uniform_width") divides the continuous data into N bins with equal width.
#' The number of bins N is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of the first random variable.
#' @param y a numeric vector of the second random variable.
#'
#' @return discretize2d.uniform_width returns a 2-dimensional count table.
#' @export
#'
#' @examples
#'
#' # two numeric vectors corresponding to two continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#'
#' # corresponding joint count table estimated by "uniform width" algorithm
#' discretize2d.uniform_width(x,y)
discretize2d.uniform_width = function(x, y)
{
  numBins1 = floor( sqrt(length(x)) )
  numBins2 = floor( sqrt(length(y)) )
  r1 = range(x)
  r2 = range(y)

  b1 = seq(from=r1[1], to=r1[2], length.out=numBins1+1 )
  b2 = seq(from=r2[1], to=r2[2], length.out=numBins2+1 )

  XY = table( cut(x, breaks=b1, include.lowest=TRUE ),
              cut(y, breaks=b2, include.lowest=TRUE ) )

  return(XY)
}


#' Discretize three sets of continuous data into 3-dimensional bins by "uniform width" method
#'
#' discretize3d.uniform_width assigns the observations of three continuous random variables to bins according to the "uniform width" method,
#' and returns a corresponding 3-dimensional count table.
#'
#' The uniform width-based method ("uniform_width") that divides the continuous data into N bins with equal width.
#' The number of bins is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of the first random variable.
#' @param y a numeric vector of the second random variable.
#' @param z a numeric vector of the third random variable.
#'
#' @return discretize3d.uniform_width returns a 3-dimensional count table.
#' @export
#'
#' @examples
#'
#' # three numeric vectors corresponding to three continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#' z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)
#'
#' # corresponding joint count table estimated by "uniform width" algorithm
#' discretize3d.uniform_width(x,y,z)
discretize3d.uniform_width = function(x, y, z)
{
  numBins1 = floor( sqrt(length(x)) )
  numBins2 = floor( sqrt(length(y)) )
  numBins3 = floor( sqrt(length(z)) )
  r1 = range(x)
  r2 = range(y)
  r3 = range(z)

  b1 = seq(from=r1[1], to=r1[2], length.out=numBins1+1 )
  b2 = seq(from=r2[1], to=r2[2], length.out=numBins2+1 )
  b3 = seq(from=r3[1], to=r3[2], length.out=numBins3+1 )

  xyz = table( cut(x, breaks=b1, include.lowest=TRUE ),
               cut(y, breaks=b2, include.lowest=TRUE ),
               cut(z, breaks=b3, include.lowest=TRUE ))
  return(xyz)
}


#' Discretize a set of continuous data into 1-dimensional bins by uniform frequency
#'
#' discretize1d.uniform_frequency assigns the observations of a continuous random variables to bins according to the "uniform frequency" method,
#' and returns a corresponding count table.
#'
#' Uniform frequency-based method ("uniform_frequency") divides the continuous data into N bins with (approximate) equal count number.
#' The number of bins N is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of a random variable.
#'
#' @return discretize1d.uniform_frequency returns a one-dimensional count table.
#' @export
#'
#' @examples
#'
#' # a numeric vector corresponding to a continuous random variable
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#'
#' # corresponding count table estimated by "uniform frequency" algorithm
#' discretize1d.uniform_frequency(x)
discretize1d.uniform_frequency = function(x)
{
  X = table( uniform_frequency(x) )

  return(X)
}


#' Discretize two sets of continuous data into 2-dimensional bins by uniform frequency
#'
#' discretize2d.uniform_frequency assigns the observations of two continuous random variables to bins according to the "uniform frequency" method,
#' and returns a corresponding 2-dimensional count table.
#'
#' Uniform frequency-based method ("uniform_frequency") divides the continuous data into N bins with (approximate) equal count number.
#' The number of bins N is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of the first random variable.
#' @param y a numeric vector of the second random variable.
#'
#' @return discretize2d.uniform_frequency returns a 2-dimensional count table.
#' @export
#'
#' @examples
#'
#' # two numeric vectors corresponding to two continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#'
#' # corresponding joint count table estimated by "uniform frequency" algorithm
#' discretize2d.uniform_frequency(x,y)
discretize2d.uniform_frequency = function(x, y)
{
  XY = table( uniform_frequency(x),
              uniform_frequency(y) )

  return(XY)
}

#' Discretize three sets of continuous data into 3-dimensional bins by uniform frequency
#'
#' discretize3d.uniform_frequency assigns the observations of three continuous random variables to bins according to the "uniform frequency" method,
#' and returns a corresponding 3-dimensional count table.
#'
#' Uniform frequency-based method ("uniform_frequency") divides the continuous data into N bins with (approximate) equal count number.
#' The number of bins N is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of the first random variable.
#' @param y a numeric vector of the second random variable.
#' @param z a numeric vector of the third random variable.
#'
#' @return discretize3d.uniform_frequency returns a 3-dimensional count table.
#' @export
#'
#' @examples
#'
#' # three numeric vectors corresponding to three continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#' z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)
#'
#' # corresponding joint count table estimated by "uniform frequency" algorithm
#' discretize3d.uniform_frequency(x,y,z)
discretize3d.uniform_frequency = function(x, y, z)
{
  xyz = table( uniform_frequency(x),
               uniform_frequency(y),
               uniform_frequency(z) )
  return(xyz)
}


##--- the function of uniform frequency ---##
uniform_frequency = function(data){
  numBins = floor(sqrt(length(data)))
  nrepl = floor(length(data)/numBins)
  nplus = sample(1:numBins, length(data) - nrepl * numBins)
  nrep = rep(nrepl, numBins)
  nrep[nplus] = nrepl + 1
  data[order(data)] = rep(seq.int(numBins), nrep)
  return(data)
}
