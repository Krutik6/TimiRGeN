#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load data
readRDS("mm_l_wikilist.rds") -> wikilist
wikilist[1] -> WP1
readRDS("elist.rds") -> elist
#perform function
WikiMatrix(e_list =  elist, wp_list = wikilist) -> wikimatrix
#internal checks
sapply(wikilist, function(x) {
sapply(elist, function(y) sum(x %in% y))}) -> wmat
#check 1
#first 5 rows of wmat should be the same as wikimatrix
test_that("wmat is mostly the same as wikimatrix", {
  expect_equal(wmat, wikimatrix[1:5,])
})
#continue
lapply(wikilist, function(x){length(x)}) -> L
rbind(wmat, Total = unlist(L)) -> wmat2
#check 2
#wmat2 and wikimatrix should be the same
test_that("wmat should be the same as wikimatrix", {
  expect_equal(wmat2, wikimatrix)
})
#savedata
saveRDS(wikimatrix, "wikimatrix.rds", compress = "xz")
