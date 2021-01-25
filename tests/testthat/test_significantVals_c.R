#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)

#load geneslist
MAE <- MultiAssayExperiment()

metadata(MAE)[["genelist"]]<- readRDS("geneslist_c.rds")

#perform SignificantVals combined method
MAE <- significantVals(MAE, method = 'c',
                       geneList = metadata(MAE)[[1]],
                       maxVal = 0.05,
                       'adjPVal')

#internal checks
filtered_genelist_man <- lapply(metadata(MAE)[[1]],
                                function(df) df[df[[grep('adjPVal',
                                names(df), value = TRUE)]] < 0.05,])

#check1
#manual output should be a list of 5
test_that("Should be a list of 5",{
    expect_equal(length(names(filtered_genelist_man)), 5)
})

#check 2
#function and manual output should be the same
test_that("manual and functional sigvals should be the same", {
    expect_equal(metadata(MAE)[[2]], filtered_genelist_man)
})

metadata(MAE)[["filtered_geneslist"]] <- filtered_genelist_man

#savedata
saveRDS(metadata(MAE)[[2]], "filtered_genelist_c.rds",compress = "xz")
