#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
miR <- mm_miR 
mRNA <- mm_mRNA 

test_that("mm_miR is a dataframe", {
  expect_true(is.data.frame(miR))
})

test_that("mm_mRNA is a dataframe", {
  expect_true(is.data.frame(mRNA))
})
