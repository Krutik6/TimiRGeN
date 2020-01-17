#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(Mfuzz)
#load data
readRDS("turnpercent.rds") -> pmat
#test function
CreateClusters(method = 'c', Percent_matrix = pmat, no.clusters = 5,
Variance = 0.5)
dev.off()
#test output
#check 1
test_that("the output is in the expected form", {
expect_equal(length(colnames(ClusterData)), 5)
expect_true(is.matrix(ClusterData))
expect_true(is.list(Clusters))
expect_equal(class(Mfuzzdata)[1], "ExpressionSet")
})
#internal checks
as.data.frame(t(pmat)) -> df
#check 2
test_that("df has the correct number of columns", {
expect_equal(length(rownames(pmat)),
length(colnames(df)))
})
#continue
df$Total <- NULL
df[vapply(df, is.factor, logical(1))] <- lapply(df[vapply(df, is.factor,
logical(1))], function(x) as.numeric(as.character(x)))
round(df, 0) -> df
na.omit(df) -> df
df -> df2
Eset <- new('ExpressionSet', exprs = as.matrix(df2))
Eset_sd <- filter.std(Eset, min.std = 0.5)
dev.off()
Eset_st <- standardise(Eset_sd)
m <- mestimate(Eset_st)
cl <- mfuzz(Eset_st, centers = 5, m=m)
#check 3
test_that("Manual and code output is the same", {
expect_equal(length(rownames(cl$membership)),
length(rownames(ClusterData)))
})
#save data
saveRDS(ClusterData, "ClusterData.rds", compress = "xz")
saveRDS(Clusters, "Clusters.rds", compress = "xz")
saveRDS(Mfuzzdata, "Mfuzzdata.rds", compress = "xz")
