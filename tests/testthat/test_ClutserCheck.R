#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(igraph)
library(MultiAssayExperiment)
#load data
readRDS("Clusters.rds") -> Clusters
#test function
ClusterCheck(Clusters = Clusters@metadata$Clusters, W = FALSE)
dev.off()
#visual check, looks fine
#check 1
test_that("should be 5 clusters", {
expect_equal(max(as.integer(rownames(Clusters@metadata$Clusters$centers))), 5)
})

