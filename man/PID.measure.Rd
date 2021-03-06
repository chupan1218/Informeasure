% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/information.measure.R
\name{PID.measure}
\alias{PID.measure}
\title{A comprehensive function for evaluating the partial information decomposition}
\usage{
PID.measure(
  XYZ,
  method = c("ML", "Jeffreys", "Laplace", "SG", "minimax", "shrink"),
  lambda.probs,
  unit = c("log", "log2", "log10"),
  verbose = TRUE
)
}
\arguments{
\item{XYZ}{a joint count distribution table of three random variables}

\item{method}{six probability estimation algorithms are available, "ML" is the default.}

\item{lambda.probs}{the shrinkage intensity, only called when the probability estimator is "shrink".}

\item{unit}{the base of the logarithm. The default is natural logarithm, which is "log".
For evaluating entropy in bits, it is suggested to set the unit to "log2".}

\item{verbose}{a logic variable. if verbose is true, report the shrinkage intensity.}
}
\value{
PID.measure returns a list that includes synergistic information, unique information from x, unique information from y,
redundant information and the sum of the four parts of information.
}
\description{
The PID.measure function is used to decompose two source information acting on the common target into four parts: joint information (synergy),
unique information from source x, unique information from source y and shared information (redundancy). The input of the PID.measure is the
joint count table.
}
\details{
Six probability estimation methods are available to evaluate the underlying bin probability from observed counts: \cr
method = "ML": maximum likelihood estimator, also referred to empirical probability, \cr
method = "Jeffreys": Dirichlet distribution estimator with prior a = 0.5, \cr
method = "Laplace": Dirichlet distribution estimator with prior a = 1, \cr
method = "SG": Dirichlet distribution estimator with prior a = 1/length(XY), \cr
method = "minimax": Dirichlet distribution estimator with prior a = sqrt(sum(XY))/length(XY), \cr
method = "shrink": shrinkage estimator.
}
\examples{
# three numeric vectors corresponding to three continuous random variables
x <- c(0.0, 0.2, 0.2, 0.7, 0.9, 0.9, 0.9, 0.9, 1.0)
y <- c(1.0, 2.0,  12, 8.0, 1.0, 9.0, 0.0, 3.0, 9.0)
z <- c(3.0, 7.0, 2.0,  11,  10,  10,  14, 2.0,  11)

# corresponding joint count table estimated by "uniform width" algorithm
XYZ <- discretize3D(x, y, z, "uniform_width")

# corresponding partial information decomposition
PID.measure(XYZ)

# corresponding count table estimated by "uniform frequency" algorithm
XYZ <- discretize3D(x, y, z, "uniform_frequency")

# corresponding partial information decomposition
PID.measure(XYZ)
}
\references{
Hausser, J., & Strimmer, K. (2009). Entropy Inference and the James-Stein Estimator, with Application to Nonlinear Gene Association Networks.
Journal of Machine Learning Research, 1469-1484.

Williams, P. L., & Beer, R. D. (2010). Nonnegative Decomposition of Multivariate Information. arXiv: Information Theory.

Chan, T. E., Stumpf, M. P., & Babtie, A. C. (2017). Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures.
Cell Systems, 5(3).
}
