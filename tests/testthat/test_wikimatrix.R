#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load data
MultiAssayExperiment() -> MAE
w_list[1:5] -> MAE@metadata$wlist
e_list -> MAE@metadata$elist
#perform function
WikiMatrix(e_list =  MAE@metadata$elist,
wp_list = MAE@metadata$wlist) -> MAE@ExperimentList$wikimatrix
#internal checks
sapply(MAE@metadata$wlist, function(x) {
sapply(MAE@metadata$elist, function(y) sum(x %in% y))}) -> wmat
#continue
lapply(MAE@metadata$wlist, function(x){length(x)}) -> L
rbind(wmat, Total = unlist(L)) -> wmat2
#check 2
#wmat2 and wikimatrix should be the same
test_that("wmat should be the same as wikimatrix", {
expect_equal(as.matrix(wmat2), as.matrix(MAE@ExperimentList$wikimatrix))
})
#savedata
saveRDS(MAE@ExperimentList$wikimatrix,
"wikimatrix.rds", compress = "xz")
