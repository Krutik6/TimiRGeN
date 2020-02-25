#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load data
MAE <- MultiAssayExperiment()
MAE@ExperimentList$Wmat <- readRDS("wikimatrix.rds")
#test TurnPercent
TurnPercent(wikiMatrix = MAE@ExperimentList$Wmat, rowInt = 6
) -> MAE@metadata$Pmat
#internal checks
as.matrix(MAE@ExperimentList$Wmat) -> df1
#check 1
#df1 should be a matrix
test_that("df1 is a matrix", {
expect_true(is.matrix(df1))
})
#continue
t(t(df1)/df1[6,]*100) -> X
format(round(X, 2), nsmall = 2) -> X
#check 2
#manual output is the same as functional ouput
test_that("check expected and observed output", {
expect_equal(MAE@metadata$Pmat, X)
})
#save file
#save data
saveRDS(MAE@metadata$Pmat, "Pmat.rds", compress = "xz")
