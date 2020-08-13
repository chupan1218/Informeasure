#' A comprehensive function for evaluating mutual information
#'
#' The MI.measure function is used to calculate the mutual information between two random variables from the joint count table.
#'
#' Six probability estimation methods are available to evaluate the underlying bin probability from observed counts: \cr
#' method = "ML": maximum likelihood estimator, also referred to empirical probability, \cr
#' method = "Jeffreys": Dirichlet distribution estimator with prior a = 0.5, \cr
#' method = "Laplace": Dirichlet distribution estimator with prior a = 1, \cr
#' method = "SG": Dirichlet distribution estimator with prior a = 1/length(XY), \cr
#' method = "minimax": Dirichlet distribution estimator with prior a = sqrt(sum(XY))/length(XY), \cr
#' method = "shrink": shrinkage estimator.
#' @param XY a joint count distribution table of two random variables.
#' @param method six probability estimation algorithms are available, "ML" is the default.
#' @param lambda.probs the shrinkage intensity, only called when the probability estimator is "shrink".
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#' @param verbose a logic variable. if verbose is true, report the shrinkage intensity.
#'
#' @return MI.measure returns the mutual information.
#' @export
#' @import entropy
#' @importFrom entropy freqs.empirical freqs.Dirichlet freqs.shrink
#'
#' @references
#'
#' Hausser, J., & Strimmer, K. (2009). Entropy Inference and the James-Stein Estimator, with Application to Nonlinear Gene Association Networks.
#' Journal of Machine Learning Research, 1469-1484.
#'
#' Wyner, A. D. (1978). A definition of conditional mutual information for arbitrary ensembles. Information & Computation, 38(1), 51-59.
#'
#' @examples
#'
#' # load Informeasure library
#' library("Informeasure")
#'
#' # two numeric vectors corresponding to two continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' XY <- discretize2d(x, y, "uniform_width")
#'
#' # corresponding mutual information
#' MI.measure(XY)
MI.measure = function(XY, method = c("ML", "Jeffreys", "Laplace", "SG", "minimax", "shrink"),
                      lambda.probs, unit = c("log", "log2", "log10"), verbose = TRUE){

  method = match.arg(method)
    unit = match.arg(unit)

  if(method == "ML"){
    probs = freqs.empirical(XY)
    MI = MI.plugin(probs, unit)
  }

  if(method == "Jeffreys"){
    probs = freqs.Dirichlet(XY, a = 1/2)
    MI = MI.plugin(probs, unit)
  }

  if(method == "Laplace"){
    probs = freqs.Dirichlet(XY, a = 1)
    MI = MI.plugin(probs, unit)
  }

  if(method == "SG"){
    probs = freqs.Dirichlet(XY, a = 1/length(XY))
    MI = MI.plugin(probs, unit)
  }

  if(method == "minimax"){
    probs = freqs.Dirichlet(XY, a = sqrt(sum(XY))/length(XY))
    MI = MI.plugin(probs, unit)
  }

  if(method == "shrink"){
    probs = freqs.shrink(XY, lambda.freqs = lambda.probs, verbose = verbose)
    MI = MI.plugin(probs, unit)
  }

  return(MI)
}


#' A comprehensive function for estimating conditional mutual information
#'
#' The CMI.measure function is used to calculate the expected mutual information between two random variables conditioned on the third one
#' from the joint count table.
#'
#' Six probability estimation methods are available to evaluate the underlying bin probability from observed counts: \cr
#' method = "ML": maximum likelihood estimator, also referred to empirical probability, \cr
#' method = "Jeffreys": Dirichlet distribution estimator with prior a = 0.5, \cr
#' method = "Laplace": Dirichlet distribution estimator with prior a = 1, \cr
#' method = "SG": Dirichlet distribution estimator with prior a = 1/length(XY), \cr
#' method = "minimax": Dirichlet distribution estimator with prior a = sqrt(sum(XY))/length(XY), \cr
#' method = "shrink": shrinkage estimator.
#' @param XYZ a joint count distribution table of three random variables.
#' @param method six probability estimation algorithms are available, "ML" is the default.
#' @param lambda.probs  the shrinkage intensity, only called when the probability estimator is "shrink".
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#' @param verbose a logic variable. if verbose is true, report the shrinkage intensity.
#'
#' @return CMI.measure returns the conditional mutual information.
#' @export
#' @import entropy
#' @importFrom entropy freqs.empirical freqs.Dirichlet freqs.shrink
#'
#' @references
#'
#' #' Hausser, J., & Strimmer, K. (2009). Entropy Inference and the James-Stein Estimator, with Application to Nonlinear Gene Association Networks.
#' Journal of Machine Learning Research, 1469-1484.
#'
#' @examples
#'
#' # load Informeasure library
#' library("Informeasure")
#'
#' # three numeric vectors corresponding to three continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#' z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' XYZ <- discretize3d(x, y, z, "uniform_width")
#'
#' # corresponding conditional mutual information
#' CMI.measure(XYZ)
CMI.measure = function(XYZ, method = c("ML", "Jeffreys", "Laplace", "SG", "minimax", "shrink"),
                       lambda.probs, unit = c("log", "log2", "log10"), verbose = TRUE){

  method = match.arg(method)
    unit = match.arg(unit)

  if(method == "ML"){
    probs = freqs.empirical(XYZ)
    CMI = CMI.plugin(probs, unit)
  }

  if(method == "Jeffreys"){
    probs = freqs.Dirichlet(XYZ, a = 1/2)
    CMI = CMI.plugin(probs, unit)
  }

  if(method == "Laplace"){
    probs = freqs.Dirichlet(XYZ, a = 1)
    CMI = CMI.plugin(probs, unit)
  }

  if(method == "SG"){
    probs = freqs.Dirichlet(XYZ, a = length(XYZ))
    CMI = CMI.plugin(probs, unit)
  }

  if(method == "minimax"){
    probs = freqs.Dirichlet(XYZ, a = sqrt(sum(XYZ))/length(XYZ))
    CMI = CMI.plugin(probs, unit)
  }

  if(method == "shrink"){
    probs = freqs.shrink(XYZ, lambda.freqs = lambda.probs, verbose = verbose)
    CMI = CMI.plugin(probs, unit)
  }

  return(CMI)
}


#' A comprehensive function for evaluating interaction information
#'
#' The II.measure function is used to calculate the amount information contained in a set of variables from the joint count table.
#' The number of variables here is limited to three.
#'
#' Six probability estimation methods are available to evaluate the underlying bin probability from observed counts: \cr
#' method = "ML": maximum likelihood estimator, also referred to empirical probability, \cr
#' method = "Jeffreys": Dirichlet distribution estimator with prior a = 0.5, \cr
#' method = "Laplace": Dirichlet distribution estimator with prior a = 1, \cr
#' method = "SG": Dirichlet distribution estimator with prior a = 1/length(XY), \cr
#' method = "minimax": Dirichlet distribution estimator with prior a = sqrt(sum(XY))/length(XY), \cr
#' method = "shrink": shrinkage estimator.
#' @param XYZ a joint count distribution table of three random variables.
#' @param method six probability estimation algorithms are available, "ML" is the default.
#' @param lambda.probs  the shrinkage intensity, only called when the probability estimator is "shrink".
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#' @param verbose a logic variable. if verbose is true, report the shrinkage intensity.
#'
#' @return II.measure returns the interaction information.
#' @export
#' @import entropy
#' @importFrom entropy freqs.empirical freqs.Dirichlet freqs.shrink
#'
#' @references
#'
#' Hausser, J., & Strimmer, K. (2009). Entropy Inference and the James-Stein Estimator, with Application to Nonlinear Gene Association Networks.
#' Journal of Machine Learning Research, 1469-1484.
#'
#' Mcgill, W. J. (1954). Multivariate information transmission. Psychometrika, 19(2), 97-116.
#'
#' @examples
#'
#' # load Informeasure library
#' library("Informeasure")
#'
#' # three numeric vectors corresponding to three continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#' z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' XYZ <- discretize3d(x, y, z, "uniform_width")
#'
#' # corresponding interaction information
#' II.measure(XYZ)
II.measure = function(XYZ, method = c("ML", "Jeffreys", "Laplace", "SG", "minimax", "shrink"),
                      lambda.probs, unit = c("log", "log2", "log10"), verbose = TRUE){

  method = match.arg(method)
    unit = match.arg(unit)

  if(method == "ML"){
    probs = freqs.empirical(XYZ)
    II = II.plugin(probs, unit)
  }

  if(method == "Jeffreys"){
    probs = freqs.Dirichlet(XYZ, a = 1/2)
    II = II.plugin(probs, unit)
  }

  if(method == "Laplace"){
    probs = freqs.Dirichlet(XYZ, a = 1)
    II = II.plugin(probs, unit)
  }

  if(method == "SG"){
    probs = freqs.Dirichlet(XYZ, a = length(XYZ))
    II = II.plugin(probs, unit)
  }

  if(method == "minimax"){
    probs = freqs.Dirichlet(XYZ, a = sqrt(sum(XYZ))/length(XYZ))
    II = II.plugin(probs, unit)
  }

  if(method == "shrink"){
    probs = freqs.shrink(XYZ, lambda.freqs = lambda.probs, verbose = verbose)
    II = II.plugin(probs, unit)
  }

  return(II)
}


#' A comprehensive function for evaluating the partial information decomposition
#'
#' The PID.measure function is used to decompose two source information acting on the common target into four parts: joint information (synergy),
#' unique information from source x, unique information from source y and shared information (redundancy). The input of the PID.measure is the
#' joint count table.
#'
#' Six probability estimation methods are available to evaluate the underlying bin probability from observed counts: \cr
#' method = "ML": maximum likelihood estimator, also referred to empirical probability, \cr
#' method = "Jeffreys": Dirichlet distribution estimator with prior a = 0.5, \cr
#' method = "Laplace": Dirichlet distribution estimator with prior a = 1, \cr
#' method = "SG": Dirichlet distribution estimator with prior a = 1/length(XY), \cr
#' method = "minimax": Dirichlet distribution estimator with prior a = sqrt(sum(XY))/length(XY), \cr
#' method = "shrink": shrinkage estimator.
#' @param XYZ a joint count distribution table of three random variables
#' @param method six probability estimation algorithms are available, "ML" is the default.
#' @param lambda.probs  the shrinkage intensity, only called when the probability estimator is "shrink".
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#' @param verbose a logic variable. if verbose is true, report the shrinkage intensity.
#'
#' @return PID.measure returns a list that includes synergistic information, unique information from x, unique information from y,
#' redundant information and the sum of the four parts of information.
#' @export
#' @import entropy
#' @importFrom entropy freqs.empirical freqs.Dirichlet freqs.shrink
#'
#' @references
#' Hausser, J., & Strimmer, K. (2009). Entropy Inference and the James-Stein Estimator, with Application to Nonlinear Gene Association Networks.
#' Journal of Machine Learning Research, 1469-1484.
#'
#' Williams, P. L., & Beer, R. D. (2010). Nonnegative Decomposition of Multivariate Information. arXiv: Information Theory.
#'
#' Chan, T. E., Stumpf, M. P., & Babtie, A. C. (2017). Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures.
#' Cell systems, 5(3).
#'
#' @examples
#'
#' # load Informeasure library
#' library("Informeasure")
#'
#' # three numeric vectors corresponding to three continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#' z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' XYZ <- discretize3d(x, y, z, "uniform_width")
#'
#' # corresponding partial information decomposition
#' PID.measure(XYZ)
#' PID.measure(XYZ, method = "Jeffreys")
#' PID.measure(XYZ, method = "Laplace")
#' PID.measure(XYZ, method = "SG")
#' PID.measure(XYZ, method = "minimax")
#' PID.measure(XYZ, method = "shrink")
#'
#' # corresponding count table estimated by "uniform frequency" algorithm
#' XYZ <- discretize3d(x, y, z, "uniform_frequency")
#'
#' # corresponding partial information decomposition
#' PID.measure(XYZ)
#' PID.measure(XYZ, method = "Jeffreys")
#' PID.measure(XYZ, method = "Laplace")
#' PID.measure(XYZ, method = "SG")
#' PID.measure(XYZ, method = "minimax")
#' PID.measure(XYZ, method = "shrink")
PID.measure = function(XYZ, method = c("ML", "Jeffreys", "Laplace", "SG", "minimax",  "shrink"),
                       lambda.probs, unit = c("log", "log2", "log10"), verbose = TRUE){

  method = match.arg(method)
    unit = match.arg(unit)

  if(method == "ML"){
    probs = freqs.empirical(XYZ)
    PID = PID.plugin(probs, unit)
  }

  if(method == "Jeffreys"){
    probs = freqs.Dirichlet(XYZ, a = 1/2)
    PID = PID.plugin(probs, unit)
  }

  if(method == "Laplace"){
    probs = freqs.Dirichlet(XYZ, a = 1)
    PID = PID.plugin(probs, unit)
  }

  if(method == "SG"){
    probs = freqs.Dirichlet(XYZ, a = length(XYZ))
    PID = PID.plugin(probs, unit)
  }

  if(method == "minimax"){
    probs = freqs.Dirichlet(XYZ, sqrt(sum(XYZ))/length(XYZ))
    PID = PID.plugin(probs, unit)
  }

  if(method == "shrink"){
    probs = freqs.shrink(XYZ, lambda.freqs = lambda.probs, verbose = verbose)
    PID = PID.plugin(probs, unit)
  }

  return(PID)
}


#' A comprehensive function for evaluating part mutual information
#'
#'
#' The PMI.measure function is used to calculate the non-linearly direct dependencies between two variables conditioned on the third one
#' form the joint count table.
#'
#' Six probability estimation methods are available to evaluate the underlying bin probability from observed counts: \cr
#' method = "ML": maximum likelihood estimator, also referred to empirical probability, \cr
#' method = "Jeffreys": Dirichlet distribution estimator with prior a = 0.5, \cr
#' method = "Laplace": Dirichlet distribution estimator with prior a = 1, \cr
#' method = "SG": Dirichlet distribution estimator with prior a = 1/length(XY), \cr
#' method = "minimax": Dirichlet distribution estimator with prior a = sqrt(sum(XY))/length(XY), \cr
#' method = "shrink": shrinkage estimator.
#' @param XYZ a joint count distribution table of three random variables.
#' @param method six probability estimation algorithms are available, "ML" is the default.
#' @param lambda.probs  the shrinkage intensity, only called when the probability estimator is "shrink".
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#' @param verbose a logic variable. if verbose is true, report the shrinkage intensity.
#'
#' @return PMI.measure returns the part mutual information.
#' @export
#' @import entropy
#' @importFrom entropy freqs.empirical freqs.Dirichlet freqs.shrink
#'
#' @references
#'
#' Hausser, J., & Strimmer, K. (2009). Entropy Inference and the James-Stein Estimator, with Application to Nonlinear Gene Association Networks.
#' Journal of Machine Learning Research, 1469-1484.
#'
#' Zhao, J., Zhou, Y., Zhang, X., & Chen, L. (2016). Part mutual information for quantifying direct associations in networks.
#' Proceedings of the National Academy of Sciences of the United States of America, 113(18), 5130-5135.
#'
#' @examples
#' # load Informeasure library
#' library("Informeasure")
#'
#' # three numeric vectors corresponding to three continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#' z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' XYZ <- discretize3d(x, y, z, "uniform_width")
#'
#' # corresponding part mutual information
#' PMI.measure(XYZ)
PMI.measure = function(XYZ, method = c("ML", "Jeffreys", "Laplace", "SG", "minimax",  "shrink"),
                       lambda.probs, unit = c("log", "log2", "log10"), verbose = TRUE){

  method = match.arg(method)
    unit = match.arg(unit)

  if(method == "ML"){
    probs = freqs.empirical(XYZ)
    PMI = PMI.plugin(probs, unit)
  }

  if(method == "Jeffreys"){
    probs = freqs.Dirichlet(XYZ, a = 1/2)
    PMI = PMI.plugin(probs, unit)
  }

  if(method == "Laplace"){
    probs = freqs.Dirichlet(XYZ, a = 1)
    PMI = PMI.plugin(probs, unit)
  }

  if(method == "SG"){
    probs = freqs.Dirichlet(XYZ, a = length(XYZ))
    PMI = PMI.plugin(probs, unit)
  }

  if(method == "minimax"){
    probs = freqs.Dirichlet(XYZ, sqrt(sum(XYZ))/length(XYZ))
    PMI = PMI.plugin(probs, unit)
  }

  if(method == "shrink"){
    probs = freqs.shrink(XYZ, lambda.freqs = lambda.probs, verbose = verbose)
    PMI = PMI.plugin(probs, unit)
  }

  return(PMI)
}
