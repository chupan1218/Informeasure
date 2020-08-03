# Informeasure
This package compiles most information measures currently available: mutual information, conditional mutual information, interaction information, partial information decomposition and part mutual information. In detail the first estimator is used to infer bivariate network while the estimation of the last four are dedicated to the identification of trivariate network. The base installation of this package allows users to approach these information measures out of the box.

# Prerequisites and Installation
```{r echo=FALSE, results='hide', message=FALSE}
library(devtools)
install_github("chupan1218/Informeasure")
```

# A simple start-up guide for using main functions of [Informeasure]
```{r echo=FALSE, results='hide', message=FALSE}
# load the package
library(Informeasure)

# input the dataset
load(system.file("extdata/tcga.brca.testdata.Rdata", package="Informeasure"))
mRNAexpression <- log2(mRNAexpression + 1)

x <- as.numeric( mRNAexpression[which(rownames(mRNAexpression)=="BRCA1"), ] )
y <- as.numeric( mRNAexpression[which(rownames(mRNAexpression)=="BARD1"), ] )

# disc
XY <- discretize2d(x,y, model = "uniform_width")

MI.measure(XY)

# 

```

# License
This project is licensed under the GPL license.

# Citation
Please cite the following paper if you use Informeasure in your research.

__*C Pan*__, YZ He, F Yang, XX Zeng* and ZL Zhang*. Informeasure: a tool to quantify the dependence between variables in biological regulatory network from an information theory perspective. **_Bio_** 2020, xx(x):x-x.

