#devtools::uses_testthat()
library(TimiRGeN)
library(MultiAssayExperiment)
library(testthat)
#load Data
readRDS("Clusters.rds") -> Clusters
#test function
ReturnCluster(ClusterData = Clusters@ExperimentList$ClusterData,
which.cluster = 1,
fit.cluster = 0.1) -> one_0_1
#check 1
test_that("returned dataframe is smaller than original", {
expect_lt(length(rownames(one_0_1)),
length(rownames(Clusters@ExperimentList$ClusterData)))
})
#continue
ReturnCluster(ClusterData = Clusters@ExperimentList$ClusterData, 
which.cluster = 1) -> one_def
ReturnCluster(ClusterData = Clusters@ExperimentList$ClusterData, 
which.cluster = 5,
fit.cluster = 0.1 ) -> five_0_1


