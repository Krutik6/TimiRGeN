#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)

#load data
MAE <- readRDS("wpdata.rds")

#run function
MAE <- reduceWiki(MAE, path_data  = assay(MAE, 3),
                  stringWiki = "Amino Acid metabolism")

#internal checks
PD <- assay(MAE, 3)
singlewiki <- PD[which(PD$name == 'Amino Acid metabolism'),]

#check 1
#test that output is the same in functional and manual code
test_that("single and singlewiki should be the same", {
    expect_equal(assay(MAE, 4), singlewiki)
})

#save data
saveRDS(singlewiki, "interactions.rds", compress = "xz")
