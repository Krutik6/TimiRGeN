setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(igraph)
#load data
readRDS("CreateClusters/Clusters.rds") -> Clusters
#test function
ClusterCheck(Clusters = Clusters, W = FALSE)
dev.off()
#visual check, looks fine
#check 1
test_that("should be 5 clusters", {
expect_equal(max(as.integer(rownames(Clusters$centers))), 5)
})

