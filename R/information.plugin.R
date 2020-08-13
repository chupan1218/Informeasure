#' A plug-in calculator for evaluating mutual information
#'
#' MI.plugin measures the mutual information between two random variables from the joint probability distribution table.
#' @param probs the joint probability distribution table of two random variables.
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#'
#' @return MI.plugin returns the mutual information.
#' @export
#' @import entropy
#' @importFrom entropy entropy.plugin
#'
#' @examples
#' # load Informeasure library
#' library("Informeasure")
#'
#' # two numeric vectors corresponding to two continuous random variables
#' x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
#' y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
#'
#' # corresponding count table estimated by "uniform width" algorithm
#' count_xy <- discretize2d(x, y, "uniform_width")
#'
#' # the joint probability distribution table of the count data
#' probs_xy <- entropy::freqs.empirical(count_xy)
#'
#' # corresponding mutual information
#' MI.plugin(probs_xy)
MI.plugin = function(probs, unit = c("log", "log2", "log10")){

  unit = match.arg(unit)

  MI = entropy.plugin(rowSums(probs),  unit = unit) +
       entropy.plugin(colSums(probs),  unit = unit) -
       entropy.plugin(probs,           unit = unit)

  #probs.x = rowSums(probs) # marginal probability
  #probs.y = colSums(probs)
  #probs.null = probs.x %o% probs.y # independence null model

  #MI = entropy::KL.plugin(probs, probs.null, unit = unit)

  return(MI)
}


#' A plug-in calculator for evaluating conditional mutual information
#'
#' CMI.plugin measures the expected mutual information between two random variables conditioned on the third one
#' from the joint probability distribution table.
#' @param probs the joint probability distribution table of three random variables.
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#'
#' @return CMI.plugin returns the conditional mutual information.
#' @export
#' @import entropy
#' @importFrom entropy entropy.plugin
#'
#' @references
#'
#' Wyner, A. D. (1978). A definition of conditional mutual information for arbitrary ensembles. Information & Computation, 38(1), 51-59.
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
#' count_xyz <- discretize3d(x, y, z, "uniform_width")
#'
#' # the joint probability distribution table of the count data
#' probs_xyz <- entropy::freqs.empirical(count_xyz)
#'
#' # corresponding conditional mutual information
#' CMI.plugin(probs_xyz)
CMI.plugin = function(probs, unit = c("log", "log2", "log10")){

  unit = match.arg(unit)

  p_XYZ = probs
  p_Z   = colSums(probs, dims = 2L)
  p_XZ  = apply(probs, c(1,3), sum)
  p_YZ  = colSums(probs, dims = 1L)

  CMI = entropy.plugin(p_XZ,  unit = unit) +
        entropy.plugin(p_YZ,  unit = unit) -
        entropy.plugin(p_XYZ, unit = unit) -
        entropy.plugin(p_Z,   unit = unit)

  return(CMI)
}


#' A plug-in calculator for evaluating the interaction information
#'
#' II.plugin measures the amount information contained in a set of variables from the joint probability distribution table.
#' The number of variables here is limited to three.
#' @param probs the joint probability distribution table of three random variables.
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#'
#' @return II.plugin returns the interaction information.
#' @export
#'
#' @references
#'
#' Mcgill, W. J. (1954). Multivariate information transmission. Psychometrika, 19(2), 97-116.
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
#' count_xyz <- discretize3d(x, y, z, "uniform_width")
#'
#' # the joint probability distribution table of the count data
#' probs_xyz <- entropy::freqs.empirical(count_xyz)
#'
#' # corresponding interaction information
#' II.plugin(probs_xyz)
II.plugin = function(probs, unit = c("log", "log2", "log10")){

  unit = match.arg(unit)

  p_XYZ = probs
  p_XY  = rowSums(probs, dims = 2L)

  II = CMI.plugin(p_XYZ, unit) - MI.plugin(p_XY, unit)

  return(II)
}


#' A plug-in calculator for evaluating partial information decomposition
#'
#' PID.plugin decomposes two source information acting on the common target into four parts: joint information (synergy),
#' unique information from source x, unique information from source y and shared information (redundancy).
#' The input of PMI.plug is the joint probability distribution table.
#' @param probs the joint probability distribution of three random variables.
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#'
#' @return PID.plugin returns a list that includes synergistic information, unique information from source x,
#' unique information from source y, redundant information and the sum of the four parts of information.
#' @export
#'
#' @references
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
#' count_xyz <- discretize3d(x, y, z, "uniform_width")
#'
#' # the joint probability distribution table of the count data
#' probs_xyz <- entropy::freqs.empirical(count_xyz)
#'
#' # corresponding partial information decomposition
#' PID.plugin(probs_xyz)
PID.plugin = function(probs, unit = c("log", "log2", "log10")){

  unit = match.arg(unit)

  p_X   = rowSums(probs, dims = 1L)
  p_XZ  = apply(probs, c(1,3), sum)
  p_Y   = apply(probs[,,], 2, sum)
  p_YZ  = colSums(probs, dims = 1L)
  p_Z   = colSums(probs, dims = 2L)

  I_XZ = MI.plugin(p_XZ, unit)
  I_YZ = MI.plugin(p_YZ, unit)

  Redundancy <- Redundancy(p_XZ, p_YZ, p_X, p_Y, p_Z, unit)

  Unique_Y = I_XZ - Redundancy
  Unique_Y = ifelse(Unique_Y>0, Unique_Y, 0)

  Unique_X = I_YZ - Redundancy
  Unique_X = ifelse(Unique_X>0, Unique_X, 0)

  Synergy = II.plugin(probs, unit) + Redundancy
  Synergy = ifelse(Synergy>0, Synergy, 0)

  PID = Synergy + Unique_Y + Unique_X + Redundancy

  return(data.frame(Synergy = Synergy,Unique_X = Unique_X, Unique_Y = Unique_Y,
                    Redundancy = Redundancy, PID = PID))
}

specific_information = function(p_iz, p_i, p_z, unit = c("log", "log2", "log10")){

  unit = match.arg(unit)

  p_i_z = t(t(p_iz) /p_z)  ##p(i|z)

  p_z_i = t(p_iz/p_i) ##p(z|i)

  tmp = t((log(1/p_z)) - (log(1/p_z_i))) * p_i_z

  if (unit == "log2")  tmp = t((log(1/p_z, 2))  - (log(1/p_z_i, 2))) * p_i_z  # change from log to log2 scale
  if (unit == "log10") tmp = t((log(1/p_z, 10)) - (log(1/p_z_i, 10))) * p_i_z # change from log to log10 scale

  specific_information = colSums(tmp, na.rm = TRUE)

  return(specific_information)
}

Redundancy = function(p_XZ, p_YZ, p_X, p_Y, p_Z, unit = c("log", "log2", "log10")){

  unit = match.arg(unit)

  specific_information_X = specific_information(p_XZ, p_X, p_Z, unit)
  specific_information_Y = specific_information(p_YZ, p_Y, p_Z, unit)
  #cat("specific_information_X: ",specific_information_X,
  #    "specific_information_y: ",specific_information_y,"\n")

  minimum_specific_information = apply(cbind(specific_information_X, specific_information_Y), 1, min)

  return(sum(p_Z * minimum_specific_information, na.rm = TRUE))
}


#' A plug-in calculator for evaluating the part mutual information
#'
#' PMI.plug measures the non-linearly direct dependencies between two variables conditioned on the third one
#' form the joint probability distribution table.
#' @param probs the joint probability distribution table of three random variables.
#' @param unit the base of the logarithm. The default is natural logarithm, which is "log".
#' For evaluating entropy in bits, it is suggested to set the unit to "log2".
#'
#' @return PMI.plugin returns the part mutual information.
#' @export
#'
#' @references
#'
#' Zhao, J., Zhou, Y., Zhang, X., & Chen, L. (2016). Part mutual information for quantifying direct associations in networks.
#' Proceedings of the National Academy of Sciences of the United States of America, 113(18), 5130-5135.
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
#' count_xyz <- discretize3d(x, y, z, "uniform_width")
#'
#' # the joint probability distribution table of the count data
#' probs_xyz <- entropy::freqs.empirical(count_xyz)
#'
#' # corresponding part mutual information
#' PMI.plugin(probs_xyz)
PMI.plugin = function(probs, unit = c("log", "log2", "log10")){

  unit = match.arg(unit)

  p_xyz = probs
  p_x   = rowSums(p_xyz, dims = 1L)
  p_y   = apply(p_xyz[,,], 2, sum)
  p_z   = colSums(p_xyz, dims = 2L)
  p_xz  = apply(p_xyz, c(1,3), sum)
  p_yz  = colSums(p_xyz, dims = 1L)

  ##--- p(x|z,y) = px_yz = p_xyz / p_yz ---##
  fun_x_yz = function(x, tmp){
    x/tmp
  }
  p_x_yz = apply(p_xyz, 1, fun_x_yz, p_yz)
  p_x_yz = array(t(p_x_yz), dim = dim(p_xyz))
  p_x_yz[is.na(p_x_yz)] = 0

  ##--- p*(x|z) = y p(x|zy)p(y) that is p_x_z = y p_x_zy * p_y ---##
  fun_x_yz_y = function(x, tmp, c){
    x*tmp[c]
  }
  p_x_z = apply(p_x_yz, 1, fun_x_yz_y, p_y, 1:length(p_y))
  p_x_z = array(t(p_x_z), dim = dim(p_xyz))
  p_x_z = array( apply(p_x_z, 3, rowSums), dim = c(dim(p_xyz)[1],1,dim(p_xyz)[3]))

  ##--- p(y|x,z) = py_xz = p_xyz / p_xz
  fun_y_xz = function(x, tmp){
    x/tmp
  }
  p_y_xz = apply(p_xyz, 2, fun_y_xz, p_xz)
  p_y_xz = array(t(p_y_xz), dim = dim(p_xyz))
  p_y_xz = array(apply(p_y_xz, 3, t), dim = dim(p_xyz))
  p_y_xz[is.na(p_y_xz)] = 0

  ##--- p*(y|z) = x p(x|zx)p(x) that is p_x_z = x p_y_xz * p_x ---##
  fun_x_yz_y = function(x, tmp, c){
    x*tmp[c]
  }
  p_y_z = apply(p_y_xz, 2, fun_x_yz_y, p_x, 1:length(p_x))
  p_y_z = array(t(p_y_z), dim = dim(p_xyz))
  p_y_z = array(apply(p_y_z, 3, t), dim = dim(p_xyz))
  p_y_z = array( colSums(p_y_z), dim = c(1, dim(p_xyz)[1], dim(p_xyz)[1]))

  ##--- p(x,y|z) = p(x,y,z)/p(z) that is p_xy_z = p_xyz/p_z ---##
  fun_xy_z = function(x, tmp, c){
    t(x)/tmp[c]
  }
  p_xy_z = apply(p_xyz, 1, fun_xy_z, p_z, 1:length(p_z))
  p_xy_z = array(t(p_xy_z), dim = dim(p_xyz))
  p_xy_z = array(t(apply(p_xy_z, 1, t)), dim = dim(p_xyz))

  ##--- p*(x|z)*p*(y|z) ----##
  tmp = as.matrix(p_x_z) %*% t(as.matrix(p_y_z))

  p_x_z_y_z = array(c(tmp[1:length(p_x_z[1,,]), 1:length(p_y_z[,,1])],
                      tmp[(1+length(p_x_z[1,,])):(2*length(p_x_z[1,,])), (1+length(p_y_z[,,1])):(2*length(p_y_z[,,1]))],
                      tmp[(1+2*length(p_x_z[1,,])):(3*length(p_x_z[1,,])), (1+2*length(p_y_z[,,1])):(3*length(p_y_z[,,1]))]), dim = dim(p_xyz))

  ##--- p(x,y|z)/(p*(x|z)*p*(y|z)) ---##
  tmp = log(p_xy_z/p_x_z_y_z)
  if (unit == "log2")   tmp = log(p_xy_z/p_x_z_y_z, 2)  # change from log to log2 scale
  if (unit == "log10")  tmp = log(p_xy_z/p_x_z_y_z, 10) # change from log to log10 scale
  tmp[is.infinite(tmp)] = 0

  PMI = sum(p_xyz * tmp, na.rm = TRUE)

  return(PMI)
}

