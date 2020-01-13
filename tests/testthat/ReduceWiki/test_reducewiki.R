setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load data
readRDS("downloadGMT/mouse_wplist.rds") -> mm_wplist
mm_wplist[[3]] -> wp2data
#run function
ReduceWiki(path_data  = wp2data,
stringWiki = "Amino Acid metabolism") -> single
#internal checks
wp2data[which(wp2data$name == 'Amino Acid metabolism')
,] -> singlewiki
#check 1
#test that output is the same in functional and manual code
test_that("single and singlewiki should be the same", {
expect_equal(single, singlewiki)
})
#save data
setwd("ReduceWiki/")
saveRDS(single, "test_net.rds", compress = "xz")
