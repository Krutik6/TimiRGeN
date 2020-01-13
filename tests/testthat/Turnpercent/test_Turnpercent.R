setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load data
readRDS("WikiMatrix/wikimatrix.rds") -> wikimatrix
#test TurnPercent
TurnPercent(wikiMatrix = wikimatrix, rowInt = 6) -> wikipercent
#internal checks
as.matrix(wikimatrix) -> df1
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
  expect_equal(wikipercent, X)
})
#save file
#save data
saveRDS(wikipercent, "Turnpercent/turnpercent.rds", compress = "xz")
