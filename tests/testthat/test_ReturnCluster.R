#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load Data
readRDS("ClusterData.rds") -> ClusterData
#test function
ReturnCluster(ClusterData = ClusterData, which.cluster = 1,
fit.cluster = 0.1) -> one_0_1
#check 1
test_that("returned dataframe is smaller than original", {
expect_lt(length(rownames(one_0_1)),
length(rownames(ClusterData)))
})
#continue
ReturnCluster(ClusterData = ClusterData, which.cluster = 1) -> one_def
ReturnCluster(ClusterData = ClusterData, which.cluster = 5,
fit.cluster = 0.1 ) -> five_0_1
#check 2
test_that("higher stringency leads to less returned networks", {
expect_lt(length(rownames(one_def)),
length(rownames(one_0_1)))
expect_equal(length(names(five_0_1)),
length(names(one_0_1)))
})

