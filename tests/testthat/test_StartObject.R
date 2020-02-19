#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
# load data
miR <- mm_miR
mRNA <- mm_mRNA
# test function
StartObject(miR = miR, mRNA = mRNA) -> MAE
Data <- list("miR" = as.data.frame(miR), "mRNA" = as.data.frame(mRNA))
MAE2 <- MultiAssayExperiment(experiments = Data)
# test 1
test_that("Functional output and manual outputs are the same",{
expect_equal(MAE, MAE2)
})
saveRDS(object = MAE, file = "MAE_mm.rds")
# human
miR <- hs_miR
mRNA <- hs_mRNA
StartObject(miR = miR, mRNA = mRNA) -> MAE
saveRDS(object = MAE, file = "MAE_hs.rds")

