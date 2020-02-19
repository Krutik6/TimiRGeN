#devtools::uses_testthat()
library(TimiRGeN)
library(MultiAssayExperiment)
library(testthat)
library(biomaRt)
#load data
mRNA <- mm_mRNA
mRNA <- mRNA[1:20,]
MAE <- StartObject(miR = NULL, mRNA = mRNA)
#test function
MAE <- getIDs_mRNA_mousetohuman(MAE = MAE, MAE@ExperimentList$mRNA,
mirror = "useast")
#check 1
test_that("mRNA output is as expected", {
expect_equal(length(names(MAE@ExperimentList$mRNA_ensembl)), 2)
expect_equal(length(names(MAE@ExperimentList$mRNA_entrez)), 2)
expect_equal(length(names(MAE@ExperimentList$mRNA_human_renamed)), 10)
})
#internal checks
musGenes <- rownames(mRNA)
human <- biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",
GRCh=37, host = "useast.ensembl.org")
mouse <- biomaRt::useEnsembl("ensembl", dataset="mmusculus_gene_ensembl",
GRCh=37, host = "useast.ensembl.org")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
values = musGenes ,
mart = mouse, attributesL = c("hgnc_symbol",
"entrezgene_id",
"ensembl_gene_id"),
martL = human, uniqueRows=TRUE)
#check 2
test_that("biomart output has four columns", {
expect_that(length(names(genesV2)), equals(4))
})
#continue
genesV2[! duplicated(genesV2$MGI.symbol),] -> mh
mh[! duplicated(mh$HGNC.symbol),] -> mhu
mhu[order(mhu$MGI.symbol),] -> mouse_human
mouse_human[which(mouse_human$MGI.symbol %in% rownames(mRNA) == TRUE),] -> mh_f
#continue
mRNA[which(rownames(mRNA) %in% mh_f$MGI.symbol == TRUE),] -> mRNA_f
mRNA_f[order(rownames(mRNA_f)),] -> mRNA_o
#check 4
test_that("mRNA_o and mh_f have the same columns",{
expect_equal(rownames(mRNA_f), mh_f$MGI.symbol)
})
#continue
cbind(mRNA_o, mh_f) -> mRNA_data
rownames(mRNA_data) <- mRNA_data$HGNC.symbol
MAE2 <- MultiAssayExperiment()
MAE2@ExperimentList$mRNA_entrez <- as.data.frame(cbind(
GENENAME = rownames(mRNA_data), ID = mRNA_data$EntrezGene.ID))
MAE2@ExperimentList$mRNA_ensembl <- as.data.frame(cbind(GENENAME =rownames(
mRNA_data), ID = mRNA_data$Gene.stable.ID))
mRNA_data$MGI.symbol <- mRNA_data$HGNC.symbol <-
mRNA_data$EntrezGene.ID <- mRNA_data$Gene.stable.ID <- NULL
MAE2@ExperimentList$mRNA_human_renamed <- mRNA_data
#check 5
test_that("manual and funcitonal outputs are the same", {
expect_equal(dim(MAE2@ExperimentList$mRNA_ensembl),
dim(MAE@ExperimentList$mRNA_ensembl))
expect_equal(dim(MAE2@ExperimentList$mRNA_entrez),
dim(MAE@ExperimentList$mRNA_entrez))
expect_equal(dim(MAE2@ExperimentList$mRNA_human_renamed),
dim(MAE@ExperimentList$mRNA_human_renamed))
})

