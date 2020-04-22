#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load data
Wmat <- readRDS("wikimatrix.rds")
MAE <- MultiAssayExperiment(list("Wikimatrix" = Wmat))
#test TurnPercent
MAE <- turnPercent(MAE, wikiMatrix = assay(MAE, 1), rowInt = 6)
#internal checks
df1 <- as.matrix(assay(MAE, 1))
#check 1
#df1 should be a matrix
test_that("df1 is a matrix", {
    expect_true(is.matrix(df1))
})
#continue
X <- t(t(df1)/df1[6,]*100)
X <- format(round(X, 2), nsmall = 2)
#check 2
#manual output is the same as functional ouput
test_that("check expected and observed output", {
    expect_equal(assay(MAE, 2), as.data.frame(X))
})
#save data
saveRDS(as.data.frame(X), "Pmat.rds", compress = "xz")
