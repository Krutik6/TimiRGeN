setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(ggplot2)
#load data
readRDS("EnrichWiki_c/EnrichWiki.rds") -> sigwiki
#use saveplots function
setwd("SavePlots/function/")
SavePlots(largeList = sigwiki, maxInt = 5, quickType = Quickdot,
          fileType = 'png')
#Internal checks
setwd("~/Documents/Package/smiRk/tests/testthat/SavePlots/test/")
plot_list = list()
for (i in 1:5) {
  p = Quickdot(X = sigwiki[[i]]@result, Y = names(sigwiki[i]))
  plot_list[[i]] = p
}
for (i in 1:5) {
  file_name = paste(names(sigwiki[i]), ".", "png", sep="")
  png(file_name)
  print(plot_list[[i]])
  dev.off()
}
setwd("~/Documents/Package/smiRk/tests/testthat/SavePlots/")
#check 1
#files made are the same
test_that("file names in test and function folder should be the same", {
  expect_equal(list.files("function/"), list.files("test/"))
})
#remove files to save space
do.call(file.remove, list(list.files("function/", full.names = TRUE)))
do.call(file.remove, list(list.files("test/", full.names = TRUE)))
#again
test_that("file names in test and function folder should be the same", {
  expect_equal(list.files("function/"), list.files("test/"))
})
