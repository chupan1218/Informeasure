testthat::test_package("Informeasure")

load(system.file("extdata/tcga.brca.testdata.Rdata", package="Informeasure"))

lncRNAexpression <- log2(lncRNAexpression + 1)
miRNAexpression  <- log2(miRNAexpression  + 1)
mRNAexpression   <- log2(mRNAexpression   + 1)

x <- as.numeric(miRNAexpression[which(rownames(miRNAexpression)   == "hsa-miR-26a-5p"), ])
y <- as.numeric(mRNAexpression[which(rownames(mRNAexpression)     == "PTEN"), ])
z <- as.numeric(lncRNAexpression[which(rownames(lncRNAexpression) == "PTENP1"), ])

XY  <- discretize2D(x,y, algorithm = "uniform_width")
XYZ <- discretize3D(x,y,z, algorithm = "uniform_frequency")

# mutual information
MI.unittesting  <- MI.measure(XY,   method = "ML",       unit = "log")

# conditional mutual information
CMI.unittesting <- CMI.measure(XYZ, method = "Jeffreys", unit = "log2")

# interaction information
II.unittesting  <- II.measure(XYZ,  method = "Laplace",  unit = "log10")

# partial information decomposition
PID.unittesting <- PID.measure(XYZ, method = "SG",       unit = "log")

# part mutual information
PMI.unittesting <- PMI.measure(XYZ, method = "minimax",  unit = "log2")
PMI.unittesting.<- PMI.measure(XYZ, method = "shrink",   unit = "log10")

test_that("unit test Informeasure", {
  expect_equal(MI.measure(XY,   method = "ML",       unit = "log"),   MI.unittesting)
  expect_equal(CMI.measure(XYZ, method = "Jeffreys", unit = "log2"),  CMI.unittesting)
  expect_equal(II.measure(XYZ,  method = "Laplace",  unit = "log10"), II.unittesting)
  expect_equal(PID.measure(XYZ, method = "SG",       unit = "log"),   PID.unittesting)
  expect_equal(PMI.measure(XYZ, method = "minimax",  unit = "log2"),  PMI.unittesting)
  expect_equal(PMI.measure(XYZ, method = "shrink",   unit = "log10"), PMI.unittesting.)
})
