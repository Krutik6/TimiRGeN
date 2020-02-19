#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(Mfuzz)
library(MultiAssayExperiment)
library(rWikiPathways)
#load data
MAE <- MultiAssayExperiment()
MAE@ExperimentList$pmat <- readRDS("Pmat.rds")
#test function
CreateClusters(MAE, method = 'c', Percent_matrix = MAE@ExperimentList$pmat,
no.clusters = 5, Variance = 0) -> MAE
dev.off()
#test output
#check 1
test_that("the output is in the expected form", {
expect_equal(length(colnames(MAE@ExperimentList$ClusterData)), 6)
expect_true(is.list(MAE@metadata$Clusters))
expect_equal(class(MAE@ExperimentList$Mfuzzdata)[1], "NULL")
})
#internal checks
as.data.frame(t(MAE@ExperimentList$pmat)) -> df
#check 2
test_that("df has the correct number of columns", {
expect_equal(length(rownames(MAE@ExperimentList$pmat)),
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
as.data.frame(cl$membership) -> X
head(X) -> X
for (i in seq_along(rownames(X))) {
getPathwayInfo(rownames(X)[i])[[3]]
} -> X$Description[i]
#check 3
test_that("Manual and code output is the same", {
expect_equal(length(colnames(X)),
length(colnames(MAE@ExperimentList$ClusterData)))
})
#save data
saveRDS(MAE, "Clusters.rds", compress = "xz")

