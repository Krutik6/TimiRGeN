#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load data
DatMat <- readRDS("MiningMatrix.rds")
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
# edit to create some interactions
DatMat$Pred_Fun[1:5] <- 1
#continue
MAE2 <- MultiAssayExperiment()
MAE2 <- matrixFilter(MAE2, miningMatrix = DatMat, predictedOnly = FALSE,
                         threshold = 1,
                         negativeOnly = FALSE)
#save data
saveRDS(assay(MAE2, 1), "filt_df.rds", compress = "xz")
