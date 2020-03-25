#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
miR <- hs_miR
mRNA <- hs_mRNA

test_that("hs_miR is a dataframe", {
  expect_true(is.data.frame(miR))
})

test_that("hs_mRNA is a dataframe", {
  expect_true(is.data.frame(mRNA))
})


