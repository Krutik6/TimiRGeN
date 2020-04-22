#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load data
MAE <- MultiAssayExperiment()
metadata(MAE)[["w_list"]] <- w_list[1:5]
metadata(MAE)[["e_list"]] <- e_list
#perform function
MAE <- wikiMatrix(MAE, ID_list =  metadata(MAE)[[2]],
                         wp_list = metadata(MAE)[[1]])
#internal checks
wmat <- sapply(metadata(MAE)[[1]], function(x) {
            sapply(metadata(MAE)[[2]],
                    function(y) sum(x %in% y))})
#continue
L <- lapply(metadata(MAE)[[1]], function(x){length(x)})
wmat2 <- rbind(wmat, Total = unlist(L))
#check 2
#wmat2 and wikimatrix should be the same
test_that("wmat should be the same as wikimatrix", {
    expect_equal(as.matrix(wmat2), as.matrix(assay(MAE, 1)))
})
#savedata
saveRDS(assay(MAE, 1),"wikimatrix.rds", compress = "xz")
