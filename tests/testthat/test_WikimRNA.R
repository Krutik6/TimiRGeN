library(biomaRt)
library(smiRk)
library(testthat)
#load data
#negative test
mm_mRNA -> mRNA
mRNA[1:200,] -> mRNA
getIDs_mRNA_mouse(mRNA)
Express(df = mRNA, dataType = 'Log2FC', genes_ID = mRNA_entrez,
idColumn = 'GENENAME') -> mRNA_express
path_data <- data.frame("wpid" = c(rep("WP571", 6)),
"gene" = c(16175, 12370,26419, 19249, 19645, 18479),
"name" = c(rep("Fas pathway and Stress induction of HSP regulation", 6)))
ReduceWiki(path_data = path_data,
stringWiki = 'Fas pathway and Stress induction of
HSP regulation') -> singlewiki
#WikimRNA(mRNA_express = mRNA_express, SingleWiki = singlewiki)
mRNA_express[which(mRNA_express$ID %in% singlewiki$gene),] -> GenesofInterest
test_that("There will be no genes of interest",{
expect_gt(1, length(rownames(GenesofInterest)))
expect_equal(length(names(GenesofInterest)),
length(names(mRNA_express)))
})
#positive test
readRDS("test_net.rds") -> test_net
readRDS("mRNA_log2fc.rds") -> mRNA_log2fc
mRNA_log2fc[which(mRNA_log2fc$ID %in% test_net$gene),] -> GenesofInterest2
test_that("There will be some genes of interest",{
expect_gt(length(rownames(GenesofInterest2)),
length(rownames(GenesofInterest)))
expect_equal(length(names(GenesofInterest2)),
length(names(mRNA_express)))
})
saveRDS(GenesofInterest2, "GenesofInterest.rds", compress = "xz")
file.remove('Rplots.pdf')
