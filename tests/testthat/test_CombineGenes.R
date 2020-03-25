#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(gtools)
library(MultiAssayExperiment)
#get data sorted
MAE <- readRDS("MAE_mm.rds")
miR <- assay(MAE, 1)[1:100,]
mRNA <- assay(MAE, 2)[1:200,]
MAE2 <- MultiAssayExperiment(list("miR2" = miR, "mRNA2" = mRNA))
MAE <- c(MAE, MAE2)
#test CombineGenes function
MAE <- CombineGenes(MAE, miR_data = assay(MAE, 3), mRNA_data = assay(MAE, 4))
#internal tests
miR_data <- assay(MAE, 3)
mRNA_data <- assay(MAE, 4)
miR_order <- gtools::mixedsort(names(miR_data))
mRNA_order <- gtools::mixedsort(names(mRNA_data)) 
#check 1
#column names in miRs and mRNAs are the same
test_that("Column names are named as expected", {
    expect_equal(as.character(names(assay(MAE, 1))),
    as.character(names(assay(MAE, 2))))
})
#continue
miR_data <- assay(MAE, 3)[miR_order]
mRNA_data <- assay(MAE, 4)[mRNA_order]
genetic_data_manual <- rbind(miR_data, mRNA_data)
# external check
#check rowsums of genetic_data_man and genetic_data
test_that("function and manual outcomes are the same", {
    expect_equal(length(rownames(genetic_data_manual)),
    length(rownames(assay(MAE, 3))) + length(rownames(assay(MAE, 4))))
})
#save data
saveRDS(assay(MAE, 5), "genetic_data.rds", compress = "xz")
