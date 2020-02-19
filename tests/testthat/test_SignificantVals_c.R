#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load geneslist
MAE <- MultiAssayExperiment()
MAE@metadata$geneslist <- readRDS("geneslist_c.rds")
#perform SignificantVals combined method
MAE@metadata$filtered_genelist <- SignificantVals(method = 'c',
geneList = MAE@metadata$geneslist, maxVal = 0.05,
'adjPVal')
#internal checks
filtered_genelist_man <- lapply(MAE@metadata$geneslist , 
function(df) df[df[[grep('adjPVal',
names(df), value = TRUE)]] < 0.05,])
#check1
#manual output should be a list of 5
test_that("Should be a list of 5",{
expect_equal(length(names(filtered_genelist_man)),
5)
})
#check 2
#each list in manual output should be less than geneslist
test_that("output should be less than input in length",{
for (i in 1:5) {
expect_gt(
length(rownames(MAE@metadata$geneslist[[i]])),
length(rownames(filtered_genelist_man[[i]]))
)}
})
#check 3
#function and manual output should be the same
test_that("manual and functional sigvals should be the same", {
expect_equal(MAE@metadata$filtered_genelist, filtered_genelist_man)
})
#savedata
saveRDS(MAE@metadata$filtered_genelist, 
"filtered_genelist_c.rds",
compress = "xz")

