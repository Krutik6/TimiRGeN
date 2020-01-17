#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load data
readRDS("MiningMatrix.rds") -> DatMat
#run function
MatrixFilter(miningMatrix = DatMat, NegativeOnly = FALSE, PredictedOnly = FALSE,
THRESHOLD = 1) -> Filt_df
#check 1
#both DFs should have the same number of columns
test_that("filt_df aspects", {
expect_equal(names(DatMat), names(Filt_df))
expect_gt(length(rownames(DatMat)),
length(rownames(Filt_df)))
})
#continue
MatrixFilter(miningMatrix = DatMat, PredictedOnly = FALSE, THRESHOLD = 1,
NegativeOnly = FALSE) -> Filt_df2
#check 2
test_that("filt_df2 is more flexible than filt_df", {
expect_equal(names(DatMat), names(Filt_df))
expect_gt(length(rownames(DatMat)),
length(rownames(Filt_df)))
})
Filt_df2[which(Filt_df2$avecor <= -0.9),]
#save data
saveRDS(Filt_df2, "filt_df.rds", compress = "xz")
