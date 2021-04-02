#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)

#load data
DatMat <- readRDS("MiningMatrix.rds")
names(DatMat)[[1]] <- "corr"
MAE <- MultiAssayExperiment()

#run function
MAE <- matrixFilter(MAE, miningMatrix = DatMat, negativeOnly = FALSE,
                        predictedOnly = FALSE, threshold = 1)

#check 1
#both DFs should have the same number of columns
test_that("filt_df aspects", {
    expect_equal(names(DatMat), names(assay(MAE, 1)))
    expect_gt(length(rownames(DatMat)),
    length(rownames(assay(MAE, 1))))
})

#continue
MAE2 <- MultiAssayExperiment()

MAE2 <- matrixFilter(MAE2, miningMatrix = DatMat, predictedOnly = FALSE,
                         threshold = 0,
                         negativeOnly = FALSE, maxCor = -0.85)
#save data
saveRDS(assay(MAE2, 1), "filt_df.rds", compress = "xz")
