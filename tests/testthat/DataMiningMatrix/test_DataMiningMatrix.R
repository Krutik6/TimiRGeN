setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
#load data
readRDS("TargetScans_data/TargetScans_results.rds") -> TargetScans_res
as.data.frame(TargetScans_res) -> TargetScans_res
readRDS("miRDB_data/miRDB_resuts.rds") -> miRDB_res
as.data.frame(miRDB_res) -> miRDB_res
readRDS("miRTarBase_data/miRTarBase_results.rds") -> miRTarBase_res
as.data.frame(miRTarBase_res) -> miRTarBase_res
readRDS("interactions_df/interactions.rds") -> Interactions
#test DataMiningMatrix
DataMiningMatrix(interactionsDF = Interactions,
targetscanInt = TargetScans_res$Targetscans_Interactions,
mirdbInt = miRDB_res$miRDB_Interactions,
mirtarbaseInt = miRTarBase_res$miRTarBase_Interactions) -> X
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
rownames(Interactions) %in% TargetScans_res$Targetscans_Interactions)) -> df2
df2 %>% mutate(miRDB = as.integer(
rownames(df2) %in% miRDB_res$miRDB_Interactions)) -> df3
df3 %>% mutate(Predicted_Interactions = TargetScan + miRDB) -> df4
df4 %>% mutate(miRTarBase = as.integer(
rownames(Interactions) %in% miRTarBase_res$miRTarBase_Interactions)) -> df5
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
expect_equal(df6, X)
})
#save data
setwd("DataMiningMatrix/")
saveRDS(X, "MiningMatrix.rds", compress = "xz")
