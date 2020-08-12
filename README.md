# Informeasure
This package compiles most information measures currently available: mutual information, conditional mutual information[1], interaction information[2], partial information decomposition[3] and part mutual information[4]. Using gene expression profile data, all these estimators can be employed to quantify the nonlinear dependence between variables in biological regulatory network inference. In detail the first estimator is used to infer bivariate network while the estimation of the last four are dedicated to the identification of trivariate network.

# Prerequisites and Installation
```{r echo=FALSE, results='hide', message=FALSE}
library(devtools)
install_github("chupan1218/Informeasure")
```

# A simple start-up guide for using main functions of Informeasure
```{r echo=FALSE, results='hide', message=FALSE}
# load the package
library(Informeasure)

# load the toy test dataset 
load(system.file("extdata/tcga.brca.testdata.Rdata", package="Informeasure"))
lncRNAexpression <- log2(lncRNAexpression + 1)
miRNAexpression  <- log2(miRNAexpression  + 1)
mRNAexpression   <- log2(mRNAexpression   + 1)

x <- as.numeric(miRNAexpression[which(rownames(miRNAexpression) == "hsa-miR-26a-5p"), ])
y <- as.numeric(mRNAexpression[which(rownames(mRNAexpression) == "PTEN"), ])
z <- as.numeric(lncRNAexpression[which(rownames(lncRNAexpression) == "PTENP1"), ])

##--- two discretization models are available: "uniform_width" and "uniform_frequency" ---##
##--- six probability estimators referred to entropy package[5] are available: "ML", "Jeffreys", "Laplace", "SG", "minimax", "shrink" ---##

# mutual information
YZ <- discretize2d(y,z, model = "uniform_width")

MI.measure(YZ, method = "ML", unit = "log")

# conditional mutual information
XYZ <- discretize3d(x, y, z, model = c("uniform_frequency"))

CMI.measure(XYZ, method = "Jeffreys", unit = "log2")

# interaction information
XYZ <- discretize3d(x, y, z, model = c("uniform_width"))

II.measure(XYZ, method = "Laplace", unit = "log10")

# partial information decomposition
XYZ <- discretize3d(x, y, z, model = c("uniform_frequency"))

PID.measure(XYZ, method = "SG", unit = "log")

# part mutual information
XYZ <- discretize3d(x, y, z, model = c("uniform_width"))

PMI.measure(XYZ, method = "minimax",  unit = "log2")

XYZ <- discretize3d(x, y, z, model = c("uniform_frequency"))

PMI.measure(XYZ, method = "shrink",   unit = "log10")

```

# License
This project is licensed under the GPL license.

# References
[1] Wyner A D. A definition of conditional mutual information for arbitrary ensembles[J]. Information & Computation, 1978, 38(1): 51-59.

[2] Mcgill W J. Multivariate information transmission[J]. Psychometrika, 1954, 19(2): 97-116. 

[3] Williams P L, Beer R D. Nonnegative Decomposition of Multivariate Information[J]. arXiv: Information Theory, 2010.

[4] Zhao J, Zhou Y, Zhang X, et al. Part mutual information for quantifying direct associations in networks[J]. Proceedings of the National Academy of Sciences of the United States of America, 2016, 113(18): 5130-5135.

[5] Hausser J. and Strimmer K. Entropy inference and the James-Stein estimator, with application to nonlinear gene association networks[J]. The Journal of Machine Learning Research, 2009, 10, 1469-1484.
