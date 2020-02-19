#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(MultiAssayExperiment)
#load data
as.data.frame(readRDS("TargetScans_results.rds")) -> TargetScans_res
as.data.frame(readRDS("miRDB_resuts.rds"))-> miRDB_res
as.data.frame(readRDS("miRTarBase_results.rds")) -> miRTarBase_res
readRDS("corrmat.rds") -> Interactions
#test DataMiningMatrix
DataMiningMatrix(corrTable = Interactions, targetscan = TargetScans_res, 
mirdb = miRDB_res, mirtarbase = miRTarBase_res) -> X
#check 1
#X should have 10 columns and the last 5 should be numeric
test_that("last 5 columns are numerics", {
expect_equal(length(names(X)), 10)
for (i in 6:10) {
expect_true(is.numeric(X[,i]))
}
})
#internal checks
Interactions %>% mutate(TargetScan = as.integer(
rownames(Interactions) %in% TargetScans_res[[1]])) -> df2
rownames(df2) <- rownames(Interactions)
df2 %>% mutate(miRDB = as.integer(rownames(df2) %in% miRDB_res[[1]])
) -> df3
rownames(df3) <- rownames(Interactions)
df3 %>% mutate(Predicted_Interactions = TargetScan + miRDB) -> df4
rownames(df4) <- rownames(Interactions)
df4 %>% mutate(miRTarBase = as.integer(rownames(df4) %in% miRTarBase_res[[1]])
) -> df5
rownames(df5) <- rownames(Interactions)
df5 %>% mutate(Pred_Fun = Predicted_Interactions + miRTarBase) -> df6
rownames(df6) <- rownames(Interactions)
#check 2
#dfs should have increasing numbers of columns
test_that("there are increasing numbers of columns", {
expect_gt(length(names(df2)),
length(names(Interactions)))
expect_gt(length(names(df3)),
length(names(df2)))
expect_gt(length(names(df4)),
length(names(df3)))
expect_gt(length(names(df5)),
length(names(df4)))
expect_gt(length(names(df6)),
length(names(df5)))
expect_equal(length(names(df6)),10)
})
#check 3
#outcomes should be the same
test_that("functional and manual output is the same", {
expect_equal(rownames(df6), rownames(X))
expect_equal(colnames(df6), colnames(X))    
})
#save data
saveRDS(X, "MiningMatrix.rds", compress = "xz")
