#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
#load geneslist
MAE <- MultiAssayExperiment()
metadata(MAE)[["geneslist"]] <- readRDS("genelist_s.rds")
#test function
MAE <- significantVals(MAE, method = 's',
                        geneList = metadata(MAE)[[1]],
                        maxVal = 0.05, stringVal = "adjPVal")
#check 1
test_that("filtered is less than original", {
    metadata(MAE)[[1]] -> a
    metadata(MAE)[[2]] -> b
    expect_gt(length(a[[1]][[1]][[1]]),
              length(b[[1]][[1]][[1]]))
})
#internal checks
filted_genelist <- lapply(metadata(MAE)[[1]], function(ls){lapply(ls,
                          function(df) df[df[[grep('adjPVal',
                                                    names(df),
                                                   value = TRUE)]] < 0.05, ])})
#check if outputs are the same
#check 2
test_that("outputs are the same", {
    expect_equal(class(metadata(MAE)[[2]]), "list")
    expect_equal(length(metadata(MAE)[[2]][[1]]), 5)
    expect_identical(filted_genelist, metadata(MAE)[[2]])
})
#save data
saveRDS(metadata(MAE)[[2]], "filtered_genelist_s.rds", compress = "xz")
