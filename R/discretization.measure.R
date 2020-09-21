#' Discretize one-dimensional continuous data into bins
#'
#' The function of discretize1D is used to assign the observations of a set of continuous random variables to bins,
#' and returns a corresponding one-dimensional count table. Two of the most common discretization methods are available:
#' "uniform width" and "uniform frequency".
#'
#' Uniform width-based method ("uniform_width") divides the continuous data into N bins with equal width.
#' Uniform frequency-based method ("uniform_frequency") divides the continuous data into N bins with (approximate) equal count number.
#' By default in both methods, the number of bins N is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of the random variable x.
#' @param model two discretization algorithms are available, "uniform_width" is the default.
#'
#' @return discretize1D returns a one-dimensional count table.
#' @export
#'
#' @examples
#'
#' # load Informeasure library
#' library("Informeasure")
#'
#' # a numeric vector corresponding to a continuous random variable
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' discretize1D(x, "uniform_width")
#'
#' # corresponding count table estimated by "uniform frequency" algorithm
#' discretize1D(x, "uniform_frequency")
discretize1D <- function(x, model = c("uniform_width", "uniform_frequency")){

  model <- match.arg(model)

  if(model == "uniform_width")
    return(discretize1d.uniform_width(x))

  if(model == "uniform_frequency")
    return(discretize1d.uniform_frequency(x))
}


#' Discretize 2-dimensional continuous data into bins
#'
#' The function of discretize2D is used to assign the observations of two sets of continuous random variables to bins,
#' and returns a corresponding two-dimensional count table. Two of the most common discretization methods are available:
#' "uniform width" and "uniform frequency".
#'
#' Uniform width-based method ("uniform_width") divides the continuous data into N bins with equal width.
#' Uniform frequency-based method ("uniform_frequency") divides the continuous data into N bins with (approximate) equal count number.
#' By default in both methods, the number of bins N is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of the random variable x.
#' @param y a numeric vector of the random variable y.
#' @param model two discretization algorithms are available, "uniform_width" is the default.
#'
#' @return discretize2D returns a 2-dimensional count table.
#' @export
#'
#' @examples
#'
#' # load Informeasure library
#' library("Informeasure")
#'
#' # two numeric vectors that correspond to two continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' discretize2D(x,y, "uniform_width")
#'
#' # corresponding count table estimated by "uniform frequency" algorithm
#' discretize2D(x,y, "uniform_frequency")
discretize2D <- function(x, y, model = c("uniform_width", "uniform_frequency")){

  model <- match.arg(model)

  if(model == "uniform_width")
    return(discretize2d.uniform_width(x,y))

  if(model == "uniform_frequency")
    return(discretize2d.uniform_frequency(x,y))
}


#' Discretize 3-dimensional continuous data into bins
#'
#' The function of discretize3D is used to assign the observations of three sets of continuous random variables to bins,
#' and returns a corresponding three-dimensional count table. Two of the most common discretization methods are available:
#' "uniform width" and "uniform frequency".
#'
#' Uniform width-based method ("uniform_width") divides the continuous data into N bins with equal width.
#' Uniform frequency-based method ("uniform_frequency") divides the continuous data into N bins with (approximate) equal count number.
#' By default in both methods, the number of bins N is initialized into a round-off value according to the square root of the data size.
#' @param x a numeric vector of the random variable x.
#' @param y a numeric vector of the random variable y.
#' @param z a numeric vector of the random variable z.
#' @param model two discretization algorithms are available, "uniform_width" is the default.
#'
#' @return discretize3D returns a 3-dimensional count table.
#' @export
#'
#' @examples
#'
#' # load Informeasure library
#' library("Informeasure")
#'
#' # three vectors that correspond to three continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#' z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' discretize3D(x,y,z, "uniform_width")
#'
#' # corresponding count table estimated by "uniform frequency" algorithm
#' discretize3D(x,y,z, "uniform_frequency")
discretize3D <- function(x, y, z, model = c("uniform_width", "uniform_frequency")){

  model <- match.arg(model)

  if(model == "uniform_width")
    return(discretize3d.uniform_width(x,y,z))

  if(model == "uniform_frequency")
    return(discretize3d.uniform_frequency(x,y,z))
}

