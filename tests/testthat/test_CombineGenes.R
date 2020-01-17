#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(gtools)
#get data sorted
mm_mRNA -> mRNA
mRNA[1:200,] -> mRNA
mm_miR -> miR
miR[1:100,] -> miR
#test CombineGenes function
CombineGenes(miR_data = miR, mRNA_data = mRNA) -> genetic_data
#internal tests
gtools::mixedsort(names(miR)) -> miR_order
gtools::mixedsort(names(mRNA)) -> mRNA_order
#check 1
#column names in miRs and mRNAs are the same
test_that("Column names are named as expected", {
expect_equal(as.character(names(miR_order)),
as.character(names(mRNA_order)))
expect_equal(mRNA_order, names(genetic_data))
expect_equal(miR_order, names(genetic_data))
})
#continue
miR_data <- miR[miR_order]
mRNA_data <- mRNA[mRNA_order]
rbind(miR_data, mRNA_data) -> genetic_data_man
# external check
#check rowsums of genetic_data_man and genetic_data
test_that("function and manual outcomes are the same", {
expect_equal(length(rownames(genetic_data_man)),
length(rownames(miR)) + length(rownames(mRNA)))
expect_equal(genetic_data, genetic_data_man)
})
#save data
saveRDS(genetic_data, "genetic_data.rds", compress = "xz")
