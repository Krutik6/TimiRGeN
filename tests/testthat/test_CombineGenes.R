#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(gtools)
library(MultiAssayExperiment)
#get data sorted
readRDS("MAE_mm.rds") -> MAE
MAE@ExperimentList$mRNA[1:200,] -> MAE@ExperimentList$mRNA
MAE@ExperimentList$miR[1:100,] -> MAE@ExperimentList$miR
#test CombineGenes function
CombineGenes(miR_data = MAE@ExperimentList$miR,
mRNA_data = MAE@ExperimentList$mRNA) -> MAE@ExperimentList$genetic_data
#internal tests
miR_data <- as.data.frame(MAE@ExperimentList$miR)
mRNA_data <- as.data.frame(MAE@ExperimentList$mRNA)
miR_order <- gtools::mixedsort(names(miR_data))
mRNA_order <- gtools::mixedsort(names(mRNA_data)) 
#check 1
#column names in miRs and mRNAs are the same
test_that("Column names are named as expected", {
expect_equal(as.character(names(miR_order)),
as.character(names(mRNA_order)))
expect_equal(mRNA_order, names(MAE@ExperimentList$genetic_data))
expect_equal(miR_order, names(MAE@ExperimentList$genetic_data))
})
#continue
miR_data <- MAE@ExperimentList$miR[miR_order]
mRNA_data <- MAE@ExperimentList$mRNA[mRNA_order]
genetic_data_manual <- rbind(miR_data, mRNA_data)
# external check
#check rowsums of genetic_data_man and genetic_data
test_that("function and manual outcomes are the same", {
expect_equal(length(rownames(genetic_data_manual)),
length(rownames(MAE@ExperimentList$miR)) + 
length(rownames(MAE@ExperimentList$mRNA)))
expect_equal(MAE@ExperimentList$genetic_data, genetic_data_manual)
})
#save data
saveRDS(MAE@ExperimentList$genetic_data, "genetic_data.rds", 
compress = "xz")
