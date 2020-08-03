# Informeasure
This package compiles most information measures currently available: mutual information, conditional mutual information, interaction information, partial information decomposition and part mutual information. In detail the first estimator is used to infer bivariate network while the estimation of the last four are dedicated to the identification of trivariate network. The base installation of this package allows users to approach these information measures out of the box.

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
##--- six probability estimators referened to entropy package are available: "ML", "Jeffreys", "Laplace", "SG", "minimax", "shrink" ---##

# mutual information
YZ <- discretize2d(y,z, model = "uniform_width")

MI.measure(YZ, method = "ML", unit = "log")

# conditional mutual information
XYZ <- discretize3d(x, y, z, model = c("uniform_frequency"))

CMI.measure(XYZ, method = "Jeffreys", unit = "log2")

# Interaction information
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

# Citation
Please cite the following paper if you use Informeasure in your research.

__*C Pan*__, YZ He, F Yang, XX Zeng* and ZL Zhang*. Informeasure: a tool to quantify the dependence between variables in biological regulatory network from an information theory perspective. **_Bio_** 2020, xx(x):x-x.

