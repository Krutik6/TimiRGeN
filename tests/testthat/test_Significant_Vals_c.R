#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load geneslist
readRDS("geneslist.rds") -> geneslist
#perform SignificantVals combined method
SignificantVals(method = 'c', geneList = geneslist, maxVal = 0.05,
'adjPVal') -> filtered_genelist
#internal checks
lapply(geneslist, function(df) df[df[[grep('adjPVal', names(df),
value = TRUE)]] < 0.05,
]) -> filtered_genelist_man
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
length(rownames(geneslist[[i]])),
length(rownames(filtered_genelist_man[[i]]))
)}
})
#check 3
#function and manual output should be the same
test_that("manual and functional sigvals should be the same", {
expect_equal(filtered_genelist, filtered_genelist_man)
})
#savedata
saveRDS(filtered_genelist, "filtered_genelist.rds",
compress = "xz")

