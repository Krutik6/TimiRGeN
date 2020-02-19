#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load geneslist
MAE <- MultiAssayExperiment()
MAE@metadata$genelist <- readRDS("genelist_s.rds")
#test function
MAE@metadata$filtered_genelist <- SignificantVals(method = 's', 
geneList = MAE@metadata$genelist, maxVal = 0.05,
stringVal = "adjPVal")
#check 1
test_that("filtered is less than original", {
expect_lt(length(MAE@metadata$filtered_genelist[[1]][[1]][[1]]),
length(MAE@metadata$genelist[[1]][[1]][[1]]))
})
#internal checks
lapply(MAE@metadata$genelist, function(ls){lapply(ls, 
function(df) df[df[[grep('adjPVal',
names(df), value = TRUE)]] < 0.05, ])}) -> filt1
#check if outputs are the same
#check 2
test_that("outputs are the same", {
expect_equal(class(filt1), "list")
expect_equal(length(filt1[[1]]), 5)
expect_identical(filt1, MAE@metadata$filtered_genelist)
})
#save data
saveRDS(MAE@metadata$filtered_genelist,
"filtered_genelist_s.rds",
compress = "xz")
