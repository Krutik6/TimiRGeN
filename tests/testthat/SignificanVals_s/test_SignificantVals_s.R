setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load geneslist
readRDS("Genelist_s/genelist_s.rds") -> genelist
#test function
SignificantVals(method = 's', geneList = genelist, maxVal = 0.05,
stringVal = "adjPVal") -> filtered_genelist
#check 1
test_that("filtered is less than original", {
expect_lt(length(filtered_genelist[[1]][[1]][[1]]),
length(genelist[[1]][[1]][[1]]))
})
#internal checks
lapply(genelist, function(ls){lapply(ls, function(df) df[df[[grep('adjPVal',
names(df), value = TRUE)]] < 0.05, ])}) -> filt1
#check if outputs are the same
#check 2
test_that("outputs are the same", {
expect_equal(class(filt1), "list")
expect_equal(length(filt1[[1]]), 5)
expect_identical(filt1, filtered_genelist)
})
#save data
saveRDS(filtered_genelist, "SignificanVals_s/filtered_genelist_s",
compress = "xz")
