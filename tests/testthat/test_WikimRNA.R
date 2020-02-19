library(biomaRt)
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load data
#negative test
readRDS("IDs_mouse_mRNA.rds") -> mRNA
Express(df = mRNA@ExperimentList$mRNA, dataType = 'Log2FC',
genes_ID = mRNA@ExperimentList$mRNA_entrez,
idColumn = 'GENENAME') -> mRNA_express
path_data <- data.frame("wpid" = c(rep("WP571", 6)),
"gene" = c(16175, 12370,26419, 19249, 19645, 18479),
"name" = c(rep("Fas pathway and Stress induction of HSP regulation", 6)))
ReduceWiki(path_data = path_data,
stringWiki = "Fas pathway and Stress induction of HSP regulation") -> singlewiki
#WikimRNA(mRNA_express = mRNA_express, SingleWiki = singlewiki)
mRNA_express[which(mRNA_express$ID %in% singlewiki$gene),] -> GenesofInterest
test_that("There will be no genes of interest",{
expect_gt(1, length(rownames(GenesofInterest)))
expect_equal(length(names(GenesofInterest)),
length(names(mRNA_express)))
})
#positive test
readRDS("interactions.rds") -> test
readRDS("log2fc.rds") -> log2fc
log2fc@ExperimentList$mRNA_log2fc[
which(log2fc@ExperimentList$mRNA_log2fc$ID %in% test$gene),
] -> GenesofInterest2

test_that("There will be some genes of interest",{
expect_gt(length(rownames(GenesofInterest2)),
length(rownames(GenesofInterest)))
expect_equal(length(names(GenesofInterest2)),
length(names(mRNA_express)))
})
saveRDS(GenesofInterest2, "GenesofInterest.rds", compress = "xz")

