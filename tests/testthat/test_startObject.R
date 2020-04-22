#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
# load data
miR <- mm_miR
mRNA <- mm_mRNA
# test function
MAE <- startObject(miR = miR, mRNA = mRNA)
Data <- list("miR" = as.data.frame(miR), "mRNA" = as.data.frame(mRNA))
MAE2 <- MultiAssayExperiment(experiments = Data)
# test 1
test_that("Functional output and manual outputs are the same",{
    expect_equal(MAE, MAE2)
})
saveRDS(MAE, file= "MAE_mm.rds", compress='xz')
# human
miR <- hs_miR
mRNA <- hs_mRNA
MAE <- startObject(miR = miR, mRNA = mRNA)

