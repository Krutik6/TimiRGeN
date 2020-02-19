#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
hs_miR -> miR
hs_mRNA -> mRNA

test_that("hs_miR is a dataframe", {
  expect_true(is.data.frame(miR))
})

test_that("hs_mRNA is a dataframe", {
  expect_true(is.data.frame(mRNA))
})


