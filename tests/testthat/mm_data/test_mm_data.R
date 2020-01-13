setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
mm_miR -> miR
mm_mRNA -> mRNA

test_that("mm_miR is a dataframe", {
  expect_true(is.data.frame(miR))
})

test_that("mm_mRNA is a dataframe", {
  expect_true(is.data.frame(mRNA))
})
